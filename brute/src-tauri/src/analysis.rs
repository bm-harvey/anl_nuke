use nukers::data_set::ArchivedDataSet;
use num::Num;

use faust::event::Event;

pub struct Analysis {}

pub trait Axis<T> {
    fn min(&self) -> f64;
    fn max(&self) -> f64;
    fn bins(&self) -> usize;
    fn title(&self) -> String;
    fn eval(&self, _thing: &T) -> Vec<f64>;
}

pub struct AxisSet<T> {
    axes: Vec<Box<dyn Axis<T>>>,
}

pub struct Module<T, R: serde::Serialize> {
    name: String,
    axes: Vec<Box<dyn Axis<T>>>,
    func: fn(&T) -> Vec<R>,
}

pub trait AnlResult: serde::Serialize {}

/*
struct MyResult {
    px: f64,
    py: f64,
    pz: f64,
    z: u8,
    a: u8,
}
impl AnlResult for MyResult {}

pub struct MyComponent {}

impl AnlComponent<Event, MyResult> for MyComponent {
    fn analyze(&self, event: &Event) -> Vec<MyResult> {
        let mut result = Vec::with_capacity(event.mult());
        for p in event.iter() {
            result.push(MyResult {
                px: p.momentum().x(),
                py: p.momentum().y(),
                pz: p.momentum().z(),
                z: p.proton_num(),
                a: p.mass_num(),
            });
        }
        result
    }
}

pub trait AnlComponent<E, R: AnlResult> {
    fn analyze(&self, event: &E) -> Vec<R>;
}

pub struct Anl<'a, E: rkyv::Archive> {
    components: Vec<Box<dyn AnlComponent<E, _>>>,
    data: &'a ArchivedDataSet<E>,
}

impl<'a, E> Anl<'a, E> {
    pub fn new(data: &'a ArchivedDataSet<E>) -> Self {
        Self {
            components: Vec::new(),
            data,
        }
    }

    pub fn attatch_component<R: AnlResult>(&mut self, module: Box<dyn AnlComponent<E, R>>) {
        self.components.push(module);
    }

    pub fn analyze(&self, event: &E) -> Vec<Vec<impl AnlResult>> {
        self.components
            .iter()
            .map(|component| component.analyze(event))
            .collect()
    }
}
*/

/* every event gets analyzed by every component,
 *
 *     _ for each event up to the buffer limit                          _
 *    |     _ for each component                                   _     |
 *    |    |                                                        |    |
 *    |    |  [ AnlResult<A>, AnlResult<B>, AnlResult<C>, ... ] ... |... |
 *    |    |                                                        |    |
 *    |    |_                                                      _|    |
 *    |_                                                                _|
 *
 */
