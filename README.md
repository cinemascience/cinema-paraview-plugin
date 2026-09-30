# Cinema ParaView Plugin

The **Cinema ParaView Plugin** provides ParaView filters, extractors,
and writers for generating and processing
[Cinema](https://cinemascience.github.io/) image databases.

The plugin represents camera configurations and image metadata
explicitly in the VTK pipeline. It supports generating **color, depth,
and data images**, processing image collections, and exporting them as
Cinema databases.

A typical workflow uses **Cinema Camera Grid** to generate camera
configurations, **Cinema Imaging** or **Cinema Color Imaging** to
produce images, optional image-processing filters, and **Cinema Writer**
to create the Cinema database.

## Building

The plugin requires ParaView, CMake, and Intel Embree. Embree 3 is used
by default; Embree 4 can be selected with `USE_EMBREE3=OFF`.

``` bash
cmake \
  -DParaView_DIR=/path/to/paraview \
  -Dembree_DIR=/path/to/embree \
  /path/to/cinema-paraview-plugin

cmake --build .
```

After building, load the **CinemaExport** plugin through ParaView's
Plugin Manager.

------------------------------------------------------------------------

## Cinema Algorithm

**Cinema Algorithm** is the common base class for the Cinema VTK
filters. It provides shared pipeline functionality used by the other
modules and is not intended to be used directly from the ParaView user
interface.

------------------------------------------------------------------------

## Cinema Camera Grid

### Purpose

**Cinema Camera Grid** generates a spherical, **phi/theta-parameterized
camera grid** around an input dataset. Camera configurations are
represented as regular VTK data, allowing them to be modified and
processed by the VTK pipeline.

### Input and Output

The input dataset defines the region around which the cameras are
generated. The output is a `vtkPointSet` whose points represent camera
positions and whose point-data arrays describe parameters such as camera
direction, view-up vector, camera height, and clipping range.

### Usage

Apply **Cinema Camera Grid** to a dataset and configure the phi/theta
camera sampling. The resulting camera dataset can be passed to **Cinema
Imaging** or used as the input of a **Cinema Color Imaging** extractor.

------------------------------------------------------------------------

## Cinema Imaging

### Purpose

**Cinema Imaging** generates depth and data images directly from
geometry using Embree ray tracing. It performs the imaging operation
entirely within the VTK pipeline and does not depend on a ParaView
Render View.

### Input and Output

The filter takes a geometric dataset and a camera dataset. It produces a
`vtkMultiBlockDataSet` containing one `vtkImageData` per camera. Images
contain depth and mapped point/cell data, while camera parameters are
stored as field data.

### Usage

Connect geometry and a camera dataset, select the image resolution, and
execute the filter. The resulting depth and data images can be processed
by other Cinema filters or written with **Cinema Writer**.

------------------------------------------------------------------------

## Cinema Color Imaging

### Purpose

**Cinema Color Imaging** generates color images by capturing an existing
ParaView Render View for each camera in a camera dataset.

It is implemented as a ParaView **extractor** because it must control
and capture a Render View rather than perform a pure VTK pipeline
transformation. This also allows rendering to work correctly in
client/server configurations.

### Input and Output

The input is a `vtkPointSet` containing camera positions together with
`CameraDir` and `CameraUp`. Orthographic projection additionally uses
`CameraHeight`.

The extractor produces a Cinema database containing one captured color
image per camera. Camera metadata and pipeline provenance are stored
alongside the images.

### Usage

Create the extractor from a camera-producing pipeline object, select the
Render View to capture, and configure the projection, resolution, output
format, compression, and database name.

Both perspective and orthographic projection are supported. The final
camera clipping range computed during rendering is recorded with each
image.

------------------------------------------------------------------------

## Cinema Image Compositing

### Purpose

**Cinema Image Compositing** combines compatible Cinema images using
their image-space data.

### Input and Output

The filter takes Cinema image data and produces composited Cinema image
data.

### Usage

Connect compatible image collections and configure the compositing
operation. The result can be processed further or exported with **Cinema
Writer**.

------------------------------------------------------------------------

## Cinema Depth Image Projection

### Purpose

**Cinema Depth Image Projection** projects depth-image samples back into
spatial coordinates using the camera information stored with the image.

### Input and Output

The input consists of Cinema depth images with their corresponding
camera metadata. The output contains the reconstructed spatial data.

### Usage

Connect compatible depth images, such as those generated by **Cinema
Imaging**, to reconstruct positions from their depth values.

------------------------------------------------------------------------

## Cinema Grid Layout

### Purpose

**Cinema Grid Layout** arranges Cinema images into a grid for
visualization and further image processing.

### Input and Output

The filter takes a collection of Cinema images and produces an image
containing the arranged grid.

### Usage

Connect a Cinema image collection and configure the desired grid layout.

------------------------------------------------------------------------

## Cinema Writer

### Purpose

**Cinema Writer** writes collections of `vtkImageData` to a Cinema
database.

It supports both PyCinema HDF5 and PNG image storage and maintains the
database's `data.csv` manifest.

### Input

A `vtkImageData` or a composite dataset containing `vtkImageData`
leaves.

Each image is treated as an individual Cinema database entry. Its field
data defines the parameter metadata associated with that image.

### Output

A Cinema database directory containing:

``` text
database.cdb/
├── data.csv
├── <hash>.h5
├── <hash>.h5
└── ...
```

or, when PNG output is selected:

``` text
database.cdb/
├── data.csv
├── <hash>.png
├── <hash>.png
└── ...
```

For PNG output, images must contain an `rgb` or `rgba` point-data array.

### Usage

Connect an image or image collection to **Cinema Writer**, choose an
output directory, output format, and compression level, and execute the
writer.

Image identity is determined from the image's field data. Field-data
names, types, shapes, and values are serialized deterministically and
used to compute the image hash.

When writing into an existing database, entries with matching hashes are
replaced while new hashes are appended. The manifest is rewritten using
the union of metadata columns present in the database.

Pipeline provenance is gathered automatically by the ParaView writer
proxy and stored as additional columns in `data.csv`.

### HDF5 Format

For HDF5 output:

-   image point-data arrays are written below `/channels`,
-   image field-data arrays are written below `/meta`, and
-   image resolution is stored as `/meta/resolution`.

Compression levels from `0` to `9` are supported. Compressed datasets
use HDF5 chunked storage.

------------------------------------------------------------------------

## Pipeline Provenance

Cinema databases written by the plugin can include provenance describing
the ParaView pipeline that produced them.

The provenance system records persistent ServerManager properties from
contributing pipeline objects and recursively follows their upstream
inputs. Internal and information-only properties are excluded.

For **Cinema Writer**, provenance is collected from the writer's input
pipeline. For **Cinema Color Imaging**, it is collected from pipelines
feeding visible representations in the selected Render View together
with the camera pipeline.

Provenance is stored as additional columns in `data.csv` and is kept
separate from the image field data used to compute image identity.
