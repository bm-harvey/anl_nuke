use std::fs::File;
use std::io::BufReader;

pub trait AnalysisModule<EventType> {
    /// Required name of the module
    fn name(&self) -> String;
    /// Runs before the event loop
    fn initialize(&mut self) {}
    /// Runs once per event
    fn analyze_event(&mut self, _event: &EventType) {}
    /// Runs after the event loop
    fn finalize(&mut self) {}
    /// Runs after finalize
    fn generate_output(&mut self, _output_directory: &str) {}
}

pub struct AnalysisManager<T> {
    modules: Vec<Box<dyn AnalysisModule<T>>>,
    input_file: Option<String>,
    output_directory: Option<String>,
    max_events: usize,
}

impl<'a, EventType: serde::Deserialize<'a>> AnalysisManager<EventType> {
    pub fn new() -> Self {
        AnalysisManager::<EventType> {
            modules: Vec::new(),
            input_file: None,
            output_directory: None,
            max_events: usize::MAX,
        }
    }

    pub fn with_module(mut self, module: Box<dyn AnalysisModule<EventType>>) -> Self {
        self.modules.push(module);
        self
    }

    pub fn with_input_file(mut self, input: &str) -> Self {
        self.input_file = Some(input.into());
        self
    }

    pub fn with_output_directory(mut self, input: &str) -> Self {
        self.output_directory = Some(input.into());
        self
    }
    
    pub fn with_max_events(mut self, max_events: usize) -> Self {
        self.max_events = max_events;
        self
    }

    pub fn input_file(&self) -> Option<&str> {
        self.input_file.as_deref()
    }

    pub fn output_directory(&self) -> Option<&str> {
        self.output_directory.as_deref()
    }

    pub fn run_analysis(mut self) {
        let input_file: String = self.input_file().expect("input file not set").into();
        let output_directory: String = self.output_directory().expect("output file not set").into();

        let in_file = File::open(input_file).unwrap();
        let in_buf = BufReader::with_capacity(100_000_000 * 8, in_file);

        self.modules
            .iter_mut()
            .for_each(|module| module.initialize());

        let mut deserializer = rmp_serde::Deserializer::new(in_buf);
        let mut event_counter = 0;
        loop {
            let event_opt = EventType::deserialize(&mut deserializer);
            match event_opt {
                Err(_) => break,
                Ok(event) => {
                    event_counter += 1;
                    if event_counter > self.max_events {
                        break;
                    }
                    if event_counter % 1000 == 0 {
                        println!("event : {}", event_counter);
                    }

                    self.modules
                        .iter_mut()
                        .for_each(|module| module.analyze_event(&event));
                }
            }
        }

        self.modules.iter_mut().for_each(|module| module.finalize());

        self.modules
            .iter_mut()
            .for_each(|module| module.generate_output(&output_directory));
    }
}

