"""Published mu-test sample, upper -> lower, mu=m_p/m_e. Frequencies are references, not asserted spin-free values."""
BENCHMARKS=[
 dict(target='6.7GHz',id='C2015-p03-2',upper='A_0_5_1_+',lower='A_0_6_0_+',frequency_MHz=6668.5192,uncertainty_Hz=800.,K_mu=-42.,frequency_source='Mueller2004 Table1; BR95',mu_source='Ellingsen2011'),
 dict(target='12.2GHz',id='C2015-p06-1',upper='E_0_2_0_E',lower='E_0_3_-1_E',frequency_MHz=12178.596415,uncertainty_Hz=135.,K_mu=-32.8,frequency_source='Cozijn2026 Table2 optical combination difference',mu_source='Bagdonaite2013 TableI'),
 dict(target='44GHz',id='benchmark-44GHz',upper='A_0_7_0_+',lower='A_0_6_1_+',frequency_MHz=44069.410,uncertainty_Hz=10000.,K_mu=5.2,frequency_source='Mueller2004 Table1; T95 measured alternative; Levshakov2022 adopted calculated44069.430MHz without a quoted uncertainty',mu_source='Levshakov2022; inverse mass-ratio convention converted'),
 dict(target='48.372GHz',id='C2015-p10-1',upper='A_0_1_0_+',lower='A_0_0_0_+',frequency_MHz=48372.4558,uncertainty_Hz=700.,K_mu=-1.,frequency_source='Mueller2004 Table2; HD73',mu_source='Bagdonaite2013 TableI'),
 dict(target='48.377GHz',id='benchmark-48.377GHz',upper='E_0_1_0_E',lower='E_0_0_0_E',frequency_MHz=48376.892,uncertainty_Hz=10000.,K_mu=-1.,frequency_source='Mueller2004 Table2; T95',mu_source='Bagdonaite2013 TableI'),
 dict(target='60.531GHz',id='benchmark-60.531GHz',upper='E_0_1_0_E',lower='E_0_2_-1_E',frequency_MHz=60531.489,uncertainty_Hz=10000.,K_mu=-7.4,frequency_source='Mueller2004 Table2',mu_source='Bagdonaite2013 TableI'),
 dict(target='95GHz',id='benchmark-95GHz',upper='A_0_8_0_+',lower='A_0_7_1_+',frequency_MHz=95169.463,uncertainty_Hz=10000.,K_mu=1.9,frequency_source='Mueller2004 Table1',mu_source='Levshakov2022; inverse mass-ratio convention converted'),
 dict(target='96.739GHz',id='benchmark-96.739GHz',upper='E_0_2_-1_E',lower='E_0_1_-1_E',frequency_MHz=96739.362,uncertainty_Hz=5000.,K_mu=-1.,frequency_source='Mueller2004 Table2',mu_source='Dapra2017'),
 dict(target='96.741GHz',id='benchmark-96.741GHz',upper='A_0_2_0_+',lower='A_0_1_0_+',frequency_MHz=96741.375,uncertainty_Hz=5000.,K_mu=-1.,frequency_source='Mueller2004 Table2',mu_source='Dapra2017'),
 dict(target='157.271GHz',id='benchmark-157.271GHz',upper='E_0_1_0_E',lower='E_0_1_-1_E',frequency_MHz=157270.851,uncertainty_Hz=10000.,K_mu=-3.5,frequency_source='Mueller2004 Table2; T95',mu_source='Dapra2017; Bagdonaite2013'),
 # Additional diagnostic is not used as an independently selected mu benchmark.
 dict(target='24.959GHz',id='C2015-p09-1',upper='E_0_5_2_E',lower='E_0_5_1_E',frequency_MHz=24959.0789,uncertainty_Hz=400.,K_mu=None,frequency_source='Mueller2004 Table1; M85',mu_source='additional profile diagnostic'),
]
for row in BENCHMARKS:row['primary_mu_sample']=row['target']!='24.959GHz'
BENCHMARKS += [
 dict(target='96.744GHz',id='context-96.744GHz',upper='E_0_2_0_E',lower='E_0_1_0_E',frequency_MHz=96744.550,uncertainty_Hz=5000.,K_mu=-1.,frequency_source='Mueller2004 Table2',mu_source='Dapra2017 full L1498 context',primary_mu_sample=False),
 dict(target='145.097GHz',id='context-145.097GHz',upper='E_0_3_-1_E',lower='E_0_2_-1_E',frequency_MHz=145097.370,uncertainty_Hz=10000.,K_mu=-1.,frequency_source='Mueller2004 Table2; T95',mu_source='Dapra2017 full L1498 context',primary_mu_sample=False),
 dict(target='145.103GHz',id='context-145.103GHz',upper='A_0_3_0_+',lower='A_0_2_0_+',frequency_MHz=145103.152,uncertainty_Hz=10000.,K_mu=-1.,frequency_source='Mueller2004 Table2; T95',mu_source='Dapra2017 full L1498 context',primary_mu_sample=False),
]
L1498_WIDTHS={'96.739GHz':46.,'96.741GHz':49.,'96.744GHz':55.,'145.097GHz':85.,'145.103GHz':77.,'157.271GHz':79.}
L1498_OBS_SIGMA={'96.739GHz':1000.,'96.741GHz':1000.,'96.744GHz':7000.,'145.097GHz':1000.,'145.103GHz':1000.,'157.271GHz':3000.}
# Reference values actually adopted by Dapra; formal errors/covariances are NOT
# added to the independent-observation-only precision comparison below.
L1498_XU_MHZ={'96.739GHz':96739.358,'96.741GHz':96741.371,'96.744GHz':96744.545,'145.097GHz':145097.435,'145.103GHz':145103.185,'157.271GHz':157270.832}
