#transfer the pre-existing TEB dataset
python transfer.py -np 199

#genrate realizaions of CMB signal
python sims.py -np 199

#calculate inverse variance map
python cal_cinv.py

#produce lensing reconstruction results
python plot.py

#estimate the results through CIB internal delensing and CIB-lensing cross delensing
python delensing.py
