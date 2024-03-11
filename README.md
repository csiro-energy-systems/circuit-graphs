# circuit-graphs

This project is a passive circuit analysis tool using a graph representation
(the [discrete maths kind][graphs], not the chart/plot kind). At present, it is
possible to create a circuit from a set of vertices and edges, and then solve
for the current through and voltage drop across each edge, the voltage at each
node relative to a sink node, and the power drawn by each component.

Circuits can consist of passive components with arbitrary impedance, which are
represented by the edges, and (ideal) transformers, which take an edge as the
primary side and a node as the source for the secondary side. Voltage sources
are represented with a source vertex, and a sink vertex must also be chosen.

Features to come include:

* Sparse matrix backend for large circuits
* Benchmark using the [IEEE 118-bus system][ieee_118]
* *(potential)* Python integration

## Usage Example

Let us create the following circuit:
```notcode
   ┌──═════───┬────═════──┐
   │    1Ω    │    0.5Ω   │
   │+         ║           ║ 
2V ◯       1Ω ║        2Ω ║
   │-         ║           ║
   │          │           │
   └──────────┴───────────┘
```

We begin by choosing the reference vertex to be that at the bottom of the circuit
diagram. Vertices have been labelled with a `u32` index for clarity, but the 
values are unimportant:
```notcode
   0          1           2
   ┌──═════───┬────═════──┐
   │    1Ω    │    0.5Ω   │
   │+         ║           ║ 
2V ◯       1Ω ║        2Ω ║
   │-         ║           ║
   │          │           │
   └──────────┴───────────┘
              3
```

Create the vertices:
```rust
let source = VertexMetadata::new(Some(2.0), 0, VertexType::VoltageSource);

let v1 = VertexMetadata::new(None, 1, VertexType::Internal);
let v2 = VertexMetadata::new(None, 2, VertexType::Internal);

let sink = VertexMetadata::new(Some(0.0), 3, VertexType::Sink);
```
 
All vertices are marked with a tag and type. Sources and sinks must have a
voltage provided while all other vertices must be instantiated with a voltage
of `None`; the code will panic if this pattern is not followed.

Now create the edges using the tags from the vertices:
```rust
let e0 = EdgeMetadata::new(0, 1, EdgeType::new_component(1.0));
let e1 = EdgeMetadata::new(1, 2, EdgeType::new_component(2.0));
let e2 = EdgeMetadata::new(2, 3, EdgeType::new_component(0.5));
let e3 = EdgeMetadata::new(1, 3, EdgeType::new_component(1.0));
```

Each edge requires the tag of its tail and head vertices, as well as the
*conductance* (not resistance) of the resistor.

Finally, create the circuit and solve the current + voltage data
```rust
let mut circuit = Circuit::new(
    vec![source, v1, v2, sink],
    vec![e0, e1, e2, e3],
);

circuit.solve_currents_and_voltages();

let solved_currents: Vec<f64> = circuit
    .graph
    .edge_weights()
    .map(|e| e.current.unwrap())
    .collect();
let solved_voltages: Vec<Option<f64>> =
    circuit.graph.node_weights().map(|e| e.voltage).collect();
   
println!("solved_currents: {:?}", &solved_currents);
println!("solved_voltages: {:?}", &solved_voltages);

// output:
// solved_currents: [1.1666666666666672, 0.33333333333333326, 0.33333333333333326, 0.8333333333333337]
// solved_voltages: [Some(2.0), Some(0.8333333333333333), Some(0.6666666666666667), Some(0.0)]
```

[graphs]: <https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)>
[ieee_118]: <https://icseg.iti.illinois.edu/ieee-118-bus-system/>
