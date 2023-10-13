# Holistic-matching-model
Sensory perception is a holistic inference process (2023). 

## Dependencies
* Philipp Berens (2022). Circular Statistics Toolbox (Directional Statistics) (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics), MATLAB Central File Exchange. 
* We have included the function needed (circ_vmpdf.m), but would like to explicitly acknowledge the use of the package.

## Code Use
* plot\_loss\_intended\_response\_bias\_std.m - Fig.2: simulate the model and plot prior, category, likelihood & posterior, loss, bias and SD.
* plot\_loss\_intended\_response\_Tomassini.m - Fig.8: simulate noisy test vs noisy probe and plot loss.
* simulate\_HMM\_noisytest\_deGardelle.m: simulate model for noisy test w/ best fitting parameters from: De Gardelle, V., Kouider, S., & Sackur, J. (2010). An oblique illusion modulated by visibility: Non-monotonic sensory integration in orientation processing. Journal of Vision, 10(10), 6-6.
* simulate\_HMM\_noisytest\_noisyprobe\_Tomassini.m: simulate model for switching test & probe w/ best fitting parameters from: Tomassini, A., Morgan, M. J., & Solomon, J. A. (2010). Orientation uncertainty reduces perceived obliquity. Vision research, 50(5), 541-547.
* simulate\_HMM\_color.m: simulate model for undelayed & delayed color estimation w/ best fitting parameters from: Bae, G. Y., Olkkonen, M., Allred, S. R., & Flombaum, J. I. (2015). Why some colors appear more memorable than others: A model combining categories and particulars in color working memory. Journal of Experimental Psychology: General, 144(4), 744.
