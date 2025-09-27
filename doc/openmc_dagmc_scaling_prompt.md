# Prompt: Implement DAGMC geometry scaling in OpenMC (like `length_multiplier` on unstructured meshes)

## Goal

Add a user-controlled **geometry scaling factor** to DAGMC geometry in OpenMC so that `.h5m` files authored in meters (or any other unit) can be scaled to OpenMC’s **centimeter** convention at runtime. The feature must:

* Add a new **length/scale multiplier** attribute for DAGMC geometry (Python API + XML).
* Apply the scaling **inside the DAGMC implementation (C++), not in `DAGMCUniverse` logic**.
* Mirror the **existing behavior** for unstructured (MOAB/libMesh) meshes (`length_multiplier`) so the UX is consistent.
* Ship with unit/integration tests and be exercised by CI.

This is the short-term implementation in OpenMC; later, parts can be upstreamed/split out to DAGMC where appropriate.

## Primary references (read these)

* **Feature request: Scale DAGMC geometry**
  “Scale geometry when making `DAGMCUniverse`” (Issue #3200). This is the direct motivation and minimal acceptance criteria. ([GitHub][1])

* **Future alignment to keep in mind:**
  “Get OpenMC to use DAGMC to automatically get cell volume information” (Issue #3334). Design with this in mind (e.g., don’t preclude querying volumes post-scale). ([GitHub][2])

* **Existing unstructured mesh scaling implementation**
  OpenMC unstructured mesh API docs (shows `length_multiplier`). ([OpenMC Documentation][3])
  Python source property wiring: `openmc.mesh.UnstructuredMesh.length_multiplier`. ([OpenMC Documentation][4])

* **OpenMC DAGMC Python API** (for where to surface the user knob)
  `openmc.DAGMCUniverse` docs + module source. ([OpenMC Documentation][5])

* **DAGMC/OpenMC usage** (background/terminology)
  DAGMC for OpenMC user guide. ([Svalinn][6])
  General OpenMC docs. ([OpenMC Documentation][7])

## Scope & placement

* **Where the scaling should live:** in the **DAGMC class** (C++ layer that actually loads/holds the DAGMC/MOAB model), **not** in the high-level `DAGMCUniverse` Python wrapper. (Per user notes and to mirror unstructured mesh behavior that applies scaling in the mesh implementation.)

* **Analogy to copy:** In `src/mesh.cpp`, `UnstructuredMesh` reads `<length_multiplier>` from XML and MOAB path applies scaling to vertex coordinates during initialization (before bounds are determined). Replicate this pattern for DAGMC geometry. (The MOAB unstructured path already scales via `mbi_->get_adjacencies` + `get_coords`/`set_coords` loops.) Use this as the implementation archetype. *(Your reference snippet included this behavior in `MOABMesh::initialize()`.)*

## User-facing API & XML

1. **Python API** (`openmc.dagmc`):

   * Add an optional keyword to `openmc.DAGMCUniverse(...)`, e.g. `length_multiplier: float | None = None`.
   * Store as an attribute and **serialize to/from XML** (see below). Keep default `None`/`1.0` behavior (no scaling).
   * Keep backwards compatibility (existing code without the argument behaves exactly the same).

2. **XML** (`geometry.xml` DAGMC entries):

   * For the `<dagmc_universe>` element (see parsing flow in `openmc/geometry` load path), add `<length_multiplier>…</length_multiplier>` in analogy with `<dagmc_universe><filename>…</filename>…</dagmc_universe>`. See how meshes write/read `length_multiplier` in HDF5/XML for naming conventions. ([OpenMC Documentation][8])

3. **Runtime behavior**:

   * At model load, the C++ DAGMC component reads the multiplier and **scales vertex coordinates** once as the geometry is constructed (before bounding box computation / acceleration structure builds).
   * Multiplier of `1.0` or omitted → no change.

## Implementation checklist (concrete)

### A) Python layer

* File: `openmc/dagmc.py`

  * Update `class DAGMCUniverse`:

    * `__init__(..., length_multiplier: float | None = None, ...)`
    * Validate numeric type (>0 recommended but allow any positive real).
    * Persist to `self.length_multiplier`.
  * XML I/O:

    * `create_xml_subelement(...)`: add subelement `<length_multiplier>` when set.
    * `from_xml_element(...)` (or equivalent constructor): parse the optional child and set `length_multiplier`.
  * Ensure `bounded_universe()/bounding_region()` remain unaffected.

* File: `openmc/geometry.py` (XML reading/writing helpers)

  * Where `<dagmc_universe>` nodes are parsed/written, ensure the new subelement survives round-trip. ([OpenMC Documentation][8])

### B) C++ layer

* Files in `src/` that handle DAGMC geometry (look for DAGMC-specific source, e.g., `dagmc.cpp` / `dagmc.hpp` or equivalent; mirror the pattern used by `MOABMesh` in `mesh.cpp` for unstructured meshes).

  * **Add a member** to hold `length_multiplier_` (default `1.0`).
  * **Read value** from XML node (sibling to filename/node information for DAGMC universes), mirroring `UnstructuredMesh::UnstructuredMesh(pugi::xml_node)` logic. (Use `check_for_node(node, "length_multiplier")` → `std::stod(get_node_value(...))` as in mesh code.)
  * **Apply scaling** after the DAGMC/MBI is created and entities loaded but **before** computing bounds or building acceleration structures:

    * Query all vertex handles (e.g., via MOAB interface).
    * For each vertex: `get_coords` → multiply by `length_multiplier_` → `set_coords`. (Same pattern as shown in the unstructured MOAB code path you pasted.)
  * **Bounds and downstream calculations** must see the scaled coordinates; call your `determine_bounds()` after scaling.
  * **HDF5/statepoint**: If there’s any serialization of DAGMC geometry metadata, include `length_multiplier` analogously to unstructured meshes’ `to_hdf5_inner()` behavior for auditability (not strictly required if DAGMC isn’t written to statepoint—follow existing conventions).

### C) Prototype & sanity via PyMOAB (optional dev aid)

Before wiring C++: verify with a tiny PyMOAB script that grabbing all vertices and multiplying coordinates yields the expected bounding box delta. This matches your note: “first in pymoab to check that it’s working like get/set coordinates” (adjacencies → verts → scale → reset).

### D) Tests

* **Python unit tests** (e.g., `tests/unit_python/test_dagmc_universe.py`)

  * Construct a trivial DAGMC `.h5m` (tiny file in test data) and instantiate `DAGMCUniverse(length_multiplier=0.01)`; verify the **reported bounding box** (via `DAGMCUniverse.bounding_box`) scales as expected vs. `length_multiplier=1.0`. ([OpenMC Documentation][5])
  * Round-trip XML: write `geometry.xml`, re-load, and confirm the `length_multiplier` value is preserved.

* **C++ tests** (where OpenMC has DAGMC tests / or add one alongside unstructured mesh tests)

  * Load the same tiny `.h5m` with and without scaling; assert bounding box numeric equality to expected scaled values and that tracking initializations succeed.
  * If feasible, a **smoke tally** (e.g., 1–10 histories, fixed source, no physics assertion—just ensure model runs).

* **Behavior parity** tests:

  * Compare **unstructured mesh scaling** code path assertions (already present) and ensure DAGMC path mirrors the same semantics (only vertices scaled; volumes/normals derived from scaled coords implicitly).

### E) CI

* **Enable GH Actions for your fork** and ensure existing OpenMC workflows run. (Your note: “Go to openmc fork—actions—select Test CI… this just runs tests on the fork.”)
* Ensure **DAGMC is enabled** in a CI job variant (OpenMC has optional DAGMC builds). Add your new tests to the **DAGMC-enabled matrix**.
* Local dev: build OpenMC from source with DAGMC ON; run `ctest` and Python tests; then push the branch and watch CI.

### F) Post-merge follow-ups (future work, not blocking)

* Consider exposing **cell volume queries** post-scale in Python (`DAGMCUniverse.get_cell_volume(id)`), in line with #3334, possibly backed by DAGMC’s volume queries or via pydagmc when present. **Design today should not preclude** that. ([GitHub][2])
* Evaluate pushing parts of the scaling hook into **DAGMC proper** so other codes can reuse.

## Acceptance criteria

* New optional user knob available in Python:
  `openmc.DAGMCUniverse(filename="model.h5m", length_multiplier=0.01)` (meters → cm) and carried through XML.
* DAGMC geometry is **actually scaled** at load time (C++), **before** bounds/accels are computed.
* **Bounding box** and **ray tracing** respect scaled geometry; no regressions with `length_multiplier=1.0` (default).
* Unit/integration tests cover:

  * XML round-trip,
  * Bounding box scaling sanity,
  * Smoke transport with scaled DAGMC.
* CI runs DAGMC-enabled tests on your fork and passes.

## Pointers & files to look at

* **Issues to satisfy/align:** #3200 (scaling) and future #3334 (volumes). ([GitHub][1])
* **Reference behavior:** Unstructured mesh `length_multiplier` (docs + Python source). ([OpenMC Documentation][3])
* **DAGMC Python wrapper:** `openmc.DAGMCUniverse` docs + module source. ([OpenMC Documentation][5])
* **Geometry XML flow:** `openmc/geometry` parsing/writing of `<dagmc_universe>` elements. ([OpenMC Documentation][8])
* **DAGMC install/usage guides** (for test setup sanity). ([Svalinn][6])

---

If you want, I can turn this into a ready-to-open GitHub issue template or a PR description, or draft the concrete Python/C++ changes and the test files next.

[1]: https://github.com/openmc-dev/openmc/issues/3200 "Scale geometry when making DAGMCUniverse · Issue #3200 · openmc-dev/openmc · GitHub"
[2]: https://github.com/openmc-dev/openmc/issues/3334 "Get OpenMC to use DAGMC to automatically get cell volume information · Issue #3334 · openmc-dev/openmc · GitHub"
[3]: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.UnstructuredMesh.html?utm_source=chatgpt.com "openmc.UnstructuredMesh"
[4]: https://docs.openmc.org/en/stable/_modules/openmc/mesh.html?utm_source=chatgpt.com "Source code for openmc.mesh"
[5]: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.DAGMCUniverse.html?utm_source=chatgpt.com "openmc.DAGMCUniverse"
[6]: https://svalinn.github.io/DAGMC/install/openmc.html?utm_source=chatgpt.com "Installing for use with OpenMC — DAGMC - Svalinn"
[7]: https://docs.openmc.org/?utm_source=chatgpt.com "The OpenMC Monte Carlo Code — OpenMC Documentation"
[8]: https://docs.openmc.org/en/v0.15.0/_modules/openmc/geometry.html?utm_source=chatgpt.com "openmc.geometry"
