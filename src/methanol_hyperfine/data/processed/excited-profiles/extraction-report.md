# Digitized Xu and Belov spectra

The experimental inputs are digitized published figures from
[Xu et al. (2019)](https://doi.org/10.1016/j.jms.2018.12.003) and
[Belov et al. (2016)](https://doi.org/10.1063/1.4954941). They are plotted traces,
not original instrument records. The source articles and supplementary
material are identified by checksums in the [manifest](manifest.json).

The manifest records each source panel, transition, calibration, frequency
allowance and extraction route. The CSV files retain source coordinates and
normalized signal values. The 18 Xu profiles contain 6,333 samples. The Belov
triplet comparison uses the 176-point profile in
`../belov-profile-check-v1/Belov-F2.csv`; its independent calibration and
selection are documented in the adjacent manifest.

## Extraction and calibration

Xu Figure 2 supplies six vector experimental paths: Q(3) through Q(7),
K=3<-2, and R(2), K=3<-2. The line geometry is read directly from the PDF,
with the R(2) path clipped to the original axes. Frequency tick coordinates
define a linear MHz-to-kHz calibration. The approximately 0.51-point stroke
corresponds to about 2.3 kHz full width. The recorded half-stroke allowance is
an extraction allowance, not a measurement standard deviation.

Xu Figures 3 and 4 are grayscale raster images of Q(7–12), K=6<-7, and
Q(8–13), K=8<-7. Gray experimental pixels are selected separately from black
fitted curves: pixels below intensity 85 are dilated by three pixels and
excluded, and gray pixels between 90 and 200 are retained. The median
experimental ordinate is sampled every two columns. Occluded pixels remain
missing. Horizontal extraction allowances are approximately 0.5–0.8 kHz;
these correlated estimates of a plotted trace are not independent raw samples.

At least two labeled frequency ticks calibrate each raster panel; additional
ticks test the affine conversion. Maximum tick-fit residuals are below 0.4 kHz.
In the Q(12) panel of Xu Figure 3, the first boundary tick is 233155.8 MHz,
followed by the labeled 233155.9-MHz tick. Belov Figure 2 uses seven ticks
from -150 to +150 kHz. Black experimental point markers are separated from
the red fitted curve. Touching markers contribute one centroid.

Source overlays were visually checked against the original figures.
Near coincident gray and black curves, selection can suppress small residuals
or emphasize the surviving excursions; this limits inference from digitized
profiles. The signal is centered on its median and divided by its experimental
range. Absolute line intensities cannot be recovered from these normalized
figures. For Belov, `source_axis_khz` is the plotted offset and
`absolute_frequency_khz` includes the nominal 259581.350-MHz transition frequency.

## Observation model

Each trace is modeled as a frequency-modulated Lamb-dip experiment with an
exact finite-modulation second-harmonic Lorentz response, a common intrinsic
HWHM and modulation amplitude, and no beam Doppler replicas. Crossovers join
pairs of dipole-allowed lines sharing an exact upper or lower hyperfine state.
Their centers are the mean component frequencies and their relative template
strengths are geometric means of the corresponding line strengths. A bounded
nonnegative amplitude accounts for uncertain saturation. Crossovers are not
intrinsic molecular transitions in the catalogues.

The thermal background uses the modulated Doppler Gaussian at 293.15 K and
molecular mass 32.026215 u, with free signal amplitude and constant/linear
baseline. The temperature is a declared room-temperature assumption because
individual scan temperatures were unavailable. Intrinsic Lorentz HWHM bounds
are retained per spectrum and must not be identified with apparent widths in
a second-derivative display.

The objective includes the 18 Xu traces and 44 retained Xu table contrasts,
without counting the plotted multiplets twice. It also includes 55 Belov
doublet splittings. The Belov triplet profile is an independent comparison,
outside the molecular fit. Empirical edge-noise estimates, a 0.04 normalized
signal floor and a cap of 50 effective samples per trace are working weights,
not a calibrated noise covariance.

The package reproduces the fit from these retained digitized inputs. Extraction
is documented by source identities, source coordinates, calibration parameters
and selection rules; the publisher PDFs and extraction environment are not
bundled. Source-specific limitations remain applicable to numerical refits.
