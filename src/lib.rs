// /// A struct representing a passive circuit. It relies on a flow graph where
// ///  - source vertices are annotated with a voltage
// ///  - edges are annotated with a weight equal to the conductance of the
// /// component (reciprocal of resistance).
// struct Circuit {
//     // Graph doesn't exist, it's just a placeholder for now
//     graph: Graph,
// }

// impl Circuit {
//     /// Create a new circuit
//     fn new() {}

//     /// Solve for the current through and voltage across each of the circuit's
//     /// resistors.
//     fn solve() {}

//     /// Get the source nodes from the graph
//     fn get_sources() {}
// }

mod circuit_graph {
    pub struct VertexMetadata<T> {
        pub voltage: T,
    }

    impl<T> VertexMetadata<T> {
        pub fn new(voltage: T) -> VertexMetadata<T> {
            VertexMetadata { voltage }
        }
    }
}

#[cfg(test)]
mod test {
    use circuit_graph::VertexMetadata;

    use crate::circuit_graph;

    #[test]
    fn create_single_vertex() {
        let v: VertexMetadata<f64> = VertexMetadata::new(1f64);

        assert!(v.voltage == 1f64);
    }

    /// Setup a simple circuit:
    /// ----------------
    /// |              |
    /// |             | |
    /// V ^           |R|
    /// |             | |
    /// |              |
    /// ----------------
    fn create_simple_circuit() {
        todo!()
    }

    /// Test that the circuit was set up properly in that the voltage of the
    /// voltage source was set to the correct value.
    #[test]
    fn check_simple_voltage_source() {
        todo!()
    }

    /// Test that the resistance of the resistor was set to the correct value.
    #[test]
    fn check_simple_resistance() {
        todo!()
    }

    /// Test that the solved current value through the resistor is correct.
    #[test]
    fn test_simple_solved_current() {
        todo!()
    }

    /// Test that the solved voltage drop across the resistor is correct.
    #[test]
    fn test_simple_solved_voltage() {
        todo!()
    }

    /// Setup a slightly more complex circuit:
    ///     __ __
    /// ----__R__---------------
    /// |          |           |
    /// |         | |         | |
    /// V         |R|         |R|
    /// |         | |         | |
    /// |          |           |
    /// ------------------------
    fn create_complex_circuit() {
        todo!()
    }

    /// Test that the circuit's voltage source was set to the correct value.
    #[test]
    fn check_complex_voltage_source() {
        todo!()
    }

    /// Test that the resistances of the circuit's resistors were set correctly.
    #[test]
    fn check_complex_resistances() {
        todo!()
    }

    /// Test that the solved current value through each resistor is correct.
    #[test]
    fn test_complex_solved_currents() {
        todo!()
    }

    /// Test that the solved voltage drop across each resistor is correct.
    #[test]
    fn complex_solved_voltages() {
        todo!()
    }
}
