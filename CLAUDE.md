This repo is an implementation of the WRF-(Chem) https://github.com/wrf-model/wrf and WRF-SFIRE (https://github.com/openwfm/WRF-SFIRE) models in .esm files.

WRF-... models are made up of process components that are composed together into the large-scale model. Implementation will take place in three stages.

In the first stage, the original fortran codes will be instrumented to output inputs and outputs from each model component during simulations. Outputs can either be instantaneous derivatives (preferred) or integrated trajectories over a time period. The code will also be instrumented to output inputs and outputs of subassemblies of components as appropriate. A select subset of these inputs and outputs will be used as tests in .esm stub files for each component or subassembly which will be added as pull requests to https://github.com/EarthSciML/EarthSciModels.

In stage 2, The stub .esm files will be filled in with the actual physics of each process so that the tests pass, then each PR will be reviewed by a human and merged. Then, the indvidual .esm components will be imported by reference into the subassembly .esm files and coupled together so that the subassembly .esm tests pass.

In stage 3, the components that have been merged into EarthSciModels will be imported by reference and coupled into top-level EqWeather.esm, EqAtmChem.esm, and EqAtmFire.esm models and tested against simulations of the full fortran models.

All model logic should be contained in .esm files. Tests and examples should be in the "tests" and "analysis" sections of .esm files, respectively. All local tests should be conducted using the EarthSciAST rust CLI binary, a copy of which should be kept untracked in the root directory of this repo.

The .esm files should be authored compositionally, with separate components, subcomponents, and expression templates used liberally and imported by reference into other components to keep the .esm code succint, simple, and human-interpretable. Avoid repeating logic or calculations, instead factor reused pieces into their own files/components/expression templates and import them by reference. Avoid the use of scripts to mechanically generate .esm files, as they often end up producing expressios that are not properly factored or succinct.

Do not store anything large in /tmp, as it is backed by RAM rather than hard drive and can cause an OOM.
