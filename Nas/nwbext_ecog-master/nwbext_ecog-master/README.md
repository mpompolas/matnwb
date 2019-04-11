# nwbext_ecog: An NWB extension for ECoG data

Author: Ben Dichter

There are three data types, `Surface`, `CorticalSurfaces`, and `ECoGSubject`. `CorticalSurfaces` is simply a group (like a folder) to put `Surface` objects into. `Surface` holds surface mesh data (vertices and triangular faces) for sections of cortex. `ECoGSubject` is an extension of `Subject` that allows you to add the `CorticalSurfaces` object to `/general/subject`.

## Usage

### python

install:
```bash
pip install git+https://github.com/bendichter/nwbext_ecog.git
```

write:
```python
import pynwb
from nwbext_ecog import CorticalSurfaces, ECoGSubject

nwbfile = pynwb.NWBFile(...)

...

cortical_surfaces = CorticalSurfaces()
## loop me
    cortical_surfaces.create_surface(name=name, faces=faces, vertices=veritices)
##
nwbfile.subject = ECoGSubject(cortical_surfaces=cortical_surfaces)
```

read:
```python
import nwbext_ecog
from pynwb import NWBHDF5IO
io = NWBHDF5IO('path_to_file.nwb','r')
nwb = io.read()
```

### MATLAB
install:
```matlab
generateExtension('/path/to/nwbext_ecog/nwbext_ecog/ecog.namespace.yaml');
```

write:
```matlab
cortical_surfaces = types.ecog.CorticalSurfaces;

%%% loop me
    surf = types.ecog.Surface('faces', faces, 'vertices', vertices);
    cortical_surfaces.surface.set(surface_name, surf);
%%%

file.subject = types.ecog.ECoGSubject(name, cortical_surfaces);
```
