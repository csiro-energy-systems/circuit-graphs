# circuit-graphs

This project is a passive circuit analysis tool using a graph representation (the [discrete maths kind][graphs], not the chart/plot kind). At present, it is possible to create a circuit from a set of nodes and edges, where every edge represents a resistor, and then solve for the current along every edge and voltage at every node. Features to come include:

* AC voltage/current support
* Support for capacitors/inductors (generalised impedance)
* *(potential)* Python integration

## Usage Example

Let us create the following circuit:

```
----__R__-------__R__---
|   1ohm   |    0.5ohm |
|+        | |         | |
V 2   1ohm|R|     2ohm|R|
|-        | |         | |
|          |           |
------------------------
```

We begin by choosing the reference node to be that at the bottom of the circuit diagram. Also, I have below labeled the nodes with an index for clarity.

```
0          1           2
----__R__-------__R__---
|   1ohm   |    0.5ohm |
|+        | |         | |
V 2   1ohm|R|     2ohm|R|
|-        | |         | |
|          |           |
------------------------
           3
```

We can now create the nodes.

```rust
let source = VertexMetadata::new(Some(2.0), 0, VertexType::Source);

let v1 = VertexMetadata::new(None, 1, VertexType::Internal);
let v2 = VertexMetadata::new(None, 2, VertexType::Internal);

let sink = VertexMetadata::new(Some(0.0), 3, VertexType::Sink);
```

Observe that all nodes are marked with a tag and type. Sources and sinks must have a voltage provided, however no matter what is supplied for an internal vertex at instantiation, it will be initialised to `None`.

With these set we can now create the edges. These use the tags we just set on the vertices.

```rust
let e0 = EdgeMetadata::new(0, 1, 1.0);
let e1 = EdgeMetadata::new(1, 2, 2.0);
let e2 = EdgeMetadata::new(2, 3, 0.5);
let e3 = EdgeMetadata::new(1, 3, 1.0);
```
Observe that each edge requires the tag of its tail and head vertices, as well as its *conductance*, not resistance. Finally we can create the circuit:

```rust
let circuit = Circuit::new(
    vec![source, v1, v2, sink],
    vec![e0, e1, e2, e3],
)
```

[graphs]: <https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)>