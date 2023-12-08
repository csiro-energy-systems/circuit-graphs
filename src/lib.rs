
/// A struct representing a passive circuit. It relies on a flow graph where
///  - source vertices are annotated with a voltage
///  - edges are annotated with a weight equal to the conductance of the 
/// component (reciprocal of resistance).
struct Circuit {
    // Graph doesn't exist, it's just a placeholder for now
    graph: Graph,
}

impl Circuit {
    /// Create a new circuit
    fn new() {}

    /// Solve for the current through and voltage across each of the circuit's
    /// resistors.
    fn solve() {}

    /// Get the source nodes from the graph
    fn get_sources() {}
}

#[cfg(tests)]
mod tests {
    /// Setup a simple circuit:
    /// ----------------
    /// |              |
    /// |             | |
    /// V ^           |R|
    /// |             | |
    /// |              |
    /// ----------------
    fn create_simple_circuit() {}

    /// Test that the circuit was set up properly in that the voltage of the
    /// voltage source was set to the correct value.
    #[test]
    fn check_simple_voltage_source() {}

    /// Test that the resistance of the resistor was set to the correct value.
    #[test]
    fn check_simple_resistance() {}

    /// Test that the solved current value through the resistor is correct.
    #[test]
    fn test_simple_solved_current() {}

    /// Test that the solved voltage drop across the resistor is correct.
    #[test]
    fn test_simple_solved_voltage() {}

    /// Setup a slightly more complex circuit:
    ///     __ __
    /// ----__R__---------------
    /// |          |           |
    /// |         | |         | |
    /// V         |R|         |R|
    /// |         | |         | |
    /// |          |           |
    /// ------------------------
    fn create_complex_circuit() {}

    /// Test that the circuit's voltage source was set to the correct value.
    #[test]
    fn check_complex_voltage_source() {}

    /// Test that the resistances of the circuit's resistors were set correctly.
    #[test]
    fn check_complex_resistances() {}

    /// Test that the solved current value through each resistor is correct.
    #[test]
    fn test_complex_solved_currents() {}

    /// Test that the solved voltage drop across each resistor is correct.
    #[test]
    fn complex_solved_voltages() {}
}
