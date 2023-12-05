use memmap2::Mmap;
use rand::Rng;
use rkyv::validation::validators::DefaultValidator;
use rkyv::{archived_root, check_archived_root};
use rkyv::{Archive, CheckBytes, Deserialize, Serialize};
use std::marker::PhantomData;

pub type ArchivedData<'a, E> = DataCollection<'a, <DataSet<E> as Archive>::Archived, E>;

/// A DataSet is a collection of events of type E. Requires that E is `rkyv::{Archive, Serialize,
/// Deserialize}` to work properly. Generally, one `DataSet<E>` is written per file, and there can be
/// many files written for an actual data set, but that can be managed later within a
/// `DataCollection<E>`
#[derive(Archive, Deserialize, Serialize)]
#[archive(check_bytes)]
pub struct DataSet<E> {
    events: Vec<E>,
}

impl<E> DataSet<E> {
    pub fn new() -> Self {
        Self { events: Vec::new() }
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    pub fn len(&self) -> usize {
        self.events.len()
    }

    pub fn clear(&mut self) {
        self.events.clear()
    }

    pub fn add_event(&mut self, event: E) {
        self.events.push(event);
    }
}

impl<E> Default for DataSet<E> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E> DataSet<E>
where
    E: Archive,
{
    pub fn read_from_rkyv(mmap: &Mmap) -> &rkyv::Archived<Self> {
        unsafe { archived_root::<DataSet<E>>(&mmap[..]) }
    }
}

impl<'a, E> DataSet<E>
where
    E: Archive,
    <E as Archive>::Archived: CheckBytes<DefaultValidator<'a>>,
{
    pub fn validated_read_from_rkyv(mmap: &'a Mmap) -> &'a rkyv::Archived<Self> {
        check_archived_root::<DataSet<E>>(&mmap[..])
            .expect("There was a problem validating the data")
    }
}

impl<E: Archive> ArchivedDataSet<E> {
    pub fn len(&self) -> usize {
        self.events.len()
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    pub fn archived_events(&self) -> &rkyv::vec::ArchivedVec<<E as Archive>::Archived> {
        &self.events
    }
    pub fn from_map(map: &Mmap) -> &ArchivedDataSet<E> {
        unsafe { rkyv::archived_root::<DataSet<E>>(&map[..]) }
    }

    pub fn from_bytes(bytes: &[u8]) -> &ArchivedDataSet<E> {
        unsafe { rkyv::archived_root::<DataSet<E>>(bytes) }
    }
}

impl<E: Archive> ArchivedDataSet<E>
//where
//<E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    pub fn archived_event_by_idx(&self, idx: usize) -> Option<&<E as Archive>::Archived> {
        if idx >= self.len() {
            None
        } else {
            Some(&self.events[idx])
        }
    }

    pub fn archived_event_by_idx_unchecked(&self, idx: usize) -> &<E as Archive>::Archived {
        &self.events[idx]
    }
}

impl<E: Archive> ArchivedDataSet<E>
where
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    pub fn event_by_idx(&self, idx: usize) -> Option<E> {
        if idx >= self.len() {
            None
        } else {
            Some(self.events[idx].deserialize(&mut rkyv::Infallible).unwrap())
        }
    }

    pub fn random_event(&self) -> Option<E> {
        self.event_by_idx(rand::thread_rng().gen())
    }
}

pub struct DataCollection<'a, D, E> {
    data_sets: Vec<&'a D>,
    mem_maps: &'a Vec<Mmap>,
    phantom: PhantomData<E>,
}

impl<'a, E: Archive> DataCollection<'a, <DataSet<E> as Archive>::Archived, E> {
    pub fn new(
        mem_maps: &'a Vec<Mmap>,
    ) -> DataCollection<'a, <DataSet<E> as Archive>::Archived, E> {
        let mut result = Self {
            data_sets: Vec::new(),
            mem_maps,
            phantom: PhantomData,
        };

        for m in result.mem_maps.iter() {
            let ds: &ArchivedDataSet<E> = unsafe { rkyv::archived_root::<DataSet<E>>(m) };
            result.data_sets.push(ds);
        }

        result
    }

    pub fn len(&self) -> usize {
        self.data_sets.iter().map(|ds| ds.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.data_sets.is_empty()
    }

    pub fn num_sets(&self) -> usize {
        self.data_sets.len()
    }

    pub fn data_set_by_idx(&self, idx: usize) -> Option<&'a ArchivedDataSet<E>> {
        if idx >= self.num_sets() {
            None
        } else {
            Some(self.data_sets[idx])
        }
    }
}

impl<'a, E: Archive> DataCollection<'a, <DataSet<E> as Archive>::Archived, E>
where
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    pub fn event_by_idx(&self, idx: usize) -> Option<E> {
        let mut starting_idx = 0;
        for ds_idx in 0..self.num_sets() {
            let ds = &self.data_sets[ds_idx];
            let events_in_ds = ds.len();
            if events_in_ds + starting_idx > idx {
                let archived = &ds.archived_events()[idx - starting_idx];

                let event: E = archived
                    .deserialize(&mut rkyv::Infallible)
                    .expect("Issue deserializing the event from its archived form");

                return Some(event);
            } else {
                starting_idx += events_in_ds;
            }
        }
        None
    }
    pub fn random_event(&self) -> E {
        self.event_by_idx(rand::thread_rng().gen_range(0, self.len()))
            .unwrap()
    }
    pub fn archived_event_by_idx(&self, idx: usize) -> Option<&<E as Archive>::Archived> {
        let mut starting_idx = 0;
        for ds_idx in 0..self.num_sets() {
            let ds = &self.data_sets[ds_idx];
            let events_in_ds = ds.len();
            if events_in_ds + starting_idx > idx {
                let archived = &ds.archived_events()[idx - starting_idx];

                return Some(archived);
            } else {
                starting_idx += events_in_ds;
            }
        }
        None
    }

    pub fn iter(&'a self) -> DataCollectionIter<'a, E> {
        DataCollectionIter::<E>::new(self)
    }

    pub fn archived_iter(&'a self) -> ArchivedDataCollectionIter<'a, E> {
        ArchivedDataCollectionIter::<E>::new(self)
    }
}

pub struct DataCollectionIter<'a, E: Archive>
where
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    data: &'a DataCollection<'a, <DataSet<E> as Archive>::Archived, E>,
    data_set: Option<&'a ArchivedDataSet<E>>,
    data_set_idx: usize,
    idx: usize,
}

impl<'a, E: Archive> DataCollectionIter<'a, E>
where
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    pub fn new(data: &'a DataCollection<'a, <DataSet<E> as Archive>::Archived, E>) -> Self {
        Self {
            data,
            data_set: data.data_set_by_idx(0),
            data_set_idx: 0,
            idx: 0,
        }
    }
}

impl<'a, E> Iterator for DataCollectionIter<'a, E>
where
    E: Archive,
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    type Item = E;

    fn next(&mut self) -> Option<Self::Item> {
        if self.data_set_idx >= self.data.num_sets() {
            return None;
        }

        self.data_set?;

        if self.idx >= self.data_set.unwrap().len() {
            self.idx = 0;
            self.data_set_idx += 1;
            self.data_set = self.data.data_set_by_idx(self.data_set_idx);
            match self.data_set {
                None => None,
                Some(_) => self.next(),
            }
        } else {
            let idx = self.idx;
            self.idx += 1;
            let archived = self.data_set.unwrap().archived_events().get(idx).unwrap();
            Some(
                archived
                    .deserialize(&mut rkyv::Infallible)
                    .expect("Issue deserializing from archived form"),
            )
        }
    }
}

pub struct ArchivedDataCollectionIter<'a, E: Archive>
where
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    data: &'a DataCollection<'a, <DataSet<E> as Archive>::Archived, E>,
    data_set: Option<&'a ArchivedDataSet<E>>,
    data_set_idx: usize,
    idx: usize,
}

impl<'a, E: Archive> ArchivedDataCollectionIter<'a, E>
where
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    pub fn new(data: &'a DataCollection<'a, <DataSet<E> as Archive>::Archived, E>) -> Self {
        Self {
            data,
            data_set: data.data_set_by_idx(0),
            data_set_idx: 0,
            idx: 0,
        }
    }
}

impl<'a, E> Iterator for ArchivedDataCollectionIter<'a, E>
where
    E: Archive,
    <E as Archive>::Archived: Deserialize<E, rkyv::Infallible>,
{
    type Item = &'a <E as Archive>::Archived;

    fn next(&mut self) -> Option<Self::Item> {
        if self.data_set_idx >= self.data.num_sets() {
            return None;
        }

        self.data_set?;

        if self.idx >= self.data_set.unwrap().len() {
            self.idx = 0;
            self.data_set_idx += 1;
            self.data_set = self.data.data_set_by_idx(self.data_set_idx);
            match self.data_set {
                None => None,
                Some(_) => self.next(),
            }
        } else {
            let idx = self.idx;
            self.idx += 1;
            Some(self.data_set.unwrap().archived_events().get(idx).unwrap())
        }
    }
}
