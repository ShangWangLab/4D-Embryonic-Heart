# 4D Embryonic Heart

This MATLAB script is designed to synchronize 4D OCT scans of the beating heart using the TDCG method. (A citation will be provided when the associated paper has been published).

## Dependencies

Developed for MATLAB version R2022a.

Required packages:
- Image Processing Toolbox
- Signal Processing Toolbox

## I/O file layout

This script takes in RAW image files specified as uint8 XYZ intensity data. The XY plane corresponds to individual B-scans, while the Z axis is the slow scanning axis. Offsets may be specified to account for file headers. Each RAW file must contain a single input channel.

Outputs:
- <name>/ - The output directory where all the TIFF files are copied to. Each image is labeled with its phase, Z position, and channel ID.
- imageIndices\_<name>.mat - a MATLAB matrix file containing the aligned image indices corresponding to each phase and Z position.
- <name>\_log.txt - A copy of all terminal outputs, when logging is enabled (verbosity >= 1).
- <name>\_figures/ - The diagnostic figures as PNG and MATLAB .fig files when figure logging is enabled (verbosity >= 1).
- <name>.ims - The TIFFs compiled as an Imaris file when the Imaris Convert path is provided.

## Example usage

A simple example. TODO.
		TDCG_sync({'structure.RAW'}, 'aligned_output', ...
			45, ... % The anticipated heart period in number of B-scans.
			\[1 1 1000/20000\], ... % Image scan scale in microns/pixel.
			\[1024 600\], ... % B-scan \[height width\]
			1:20000 ... % The subset of B-scans to use, typically all of them.
		);

See [alignment_example.m](https://github.com/ShangWangLab/4D-Embryonic-Heart/alignment_example.m) for a more complex example.

## License

See the [LICENSE](https://github.com/ShangWangLab/4D-Embryonic-Heart/LICENSE) for more information.