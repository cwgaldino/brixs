# %% --------------------- Example 9: peak finding ------------------------ %% #
x = np.linspace(0, 10, 100)
y = br.gaussian_fwhm(x=x, amp=4, c=4, w=1) + br.gaussian_fwhm(x=x, amp=2, c=8, w=1)
s = br.Spectrum(x, y)

# uses scipy.signal.find_peaks()
s.find_peaks()
print(s.peaks)  # print list of peaks found

print(s.peaks[0])
print(s.peaks[1])

# peaks can be removed
s.peaks.remove(0)

# peaks can be added
s.peaks.append(amp=4, c=4, w=1)

# peaks can be modifies
s.peaks[1]['amp'].value = 3

# one can calculate a spectrum from peaks
s2 = s.peaks.calculate_spectrum()

# quick plot
br.figure()  
s.plot(marker='o') 
s.peaks.plot()     # place a marker at every peak found  
s2.plot()

# %% --------------------- Example 10: peak fitting ----------------------- %% #
x = np.linspace(0, 10, 100)
y = br.gaussian_fwhm(x=x, amp=4, c=4, w=1) + br.gaussian_fwhm(x=x, amp=2, c=8, w=1)
s = br.Spectrum(x, y)

# fit one peak
s.fit_peak(ranges=(0, 5))
print(s.peaks)
one_peak = s.peaks.calculate_spectrum()

# find all peaks and fit
s.find_peaks()
s.fit_peaks() 
fit = s.peaks.calculate_spectrum()

# quick plot
br.figure()  
s.plot(marker='o') 
fit.plot()
one_peak.plot()   

# %% ----------------- Example 11: defining new methods ------------------- %% #
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
s1 = br.Spectrum(x, y)

# define new function that acts on x and y
def new_function(self, value):
    """Adds a value to the y-axis"""
    self.y = self.y + value

# attaches function to type Spectrum
br.Spectrum.new_function = new_function

# from now one, every spectrum type will have new_function defined.
# even spectrum types created before defying new_function
s2 = br.Spectrum(x, y)

s1.new_function(5)
print(s1.y)

s2.new_function(2)
print(s2.y)












# %% --------------------- Example 10: peak finding ----------------------- %% #
x  = np.linspace(0, 10, 100)
y  = br.gaussian_fwhm(x=x, amp=4, c=4, w=1) + br.gaussian_fwhm(x=x, amp=2, c=8, w=1)
s1 = br.Spectrum(x, y)

x  = np.linspace(0, 10, 100)
y  = br.gaussian_fwhm(x=x, amp=5, c=2, w=1) + br.gaussian_fwhm(x=x, amp=5, c=6, w=1)
s2 = br.Spectrum(x, y)

ss = br.Spectra(s1, s2)

# uses scipy.signal.find_peaks()
ss.find_peaks()
print(ss.peaks)  # print list of peaks found
# note that peak parameters have 2 indexes, peak number and spectrum number

print(ss.peaks[0])  # peaks for spectrum 0
print(ss.peaks[1])  # peaks for spectrum 1

# quick plot
br.figure()  
lines  = ss.plot(vi=-100, marker='o')     # plot cascaded lines
offset = [line.offset for line in lines]  # get offset from cascaded lines
ss.peaks.plot(offset=offset)              # place a marker at every peak found  

# note that the ss.find_peaks() functions define peaks at a Spectra level
print(ss.peaks)

################## TODO from here

# AND at a Spectrum level
print(ss[0].peaks)
print(ss[1].peaks)

# however, the peaks from a Spectra level and Spectrum level do NOT point at the same object
# for instance, if one make changes on the ss.peaks, it does not transfer to ss[i].peaks
print(ss.peaks[0][0]['amp'].value)
print(ss[0].peaks[0]['amp'].value)
ss.peaks[0]['amp_0'].value = 444

ss[0].peaks._get_peaks_by_index(i1=0)
ss.peaks._get_peaks_by_index(i1=0, i2=0)

for a in ss.peaks:  
    print(a)

# this works! -> spectrum peaks
a = ss[0].peaks._get_peaks_by_index(i1=0)
print(a['amp_0'].value)
a['amp_0'].value = 15
print(a['amp_0'].value)
a = ss[0].peaks._get_peaks_by_index(i1=0)
print(a['amp_0'].value)

# this also works! -> spectra peaks
a = ss.peaks._get_peaks_by_index(i1=0, i2=0)
print(a['amp_0_0'].value)
a['amp_0_0'].value = 12
print(a['amp_0_0'].value)
a = ss.peaks._get_peaks_by_index(i1=0, i2=0)
print(a['amp_0_0'].value)


ss.peaks['amp_0_0'].value
ss.peaks['amp_0_0'].value = 4.2


a = {'amp_0_0': ss.peaks['amp_0_0'], }
a['amp_0_0'].value = 5.2
a['amp_0_0'].value

a = ss.peaks._get_peaks_by_index(i1=0, i2=0)


print(a['amp_0'].value)
a['amp_0'].value = 10
print(a['amp_0'].value)
print(ss.peaks[0]['amp'].value)



[0]['amp'].value = 4.2
print(ss.peaks[0][0]['amp'].value)
print(ss[0].peaks[0]['amp'].value)

# copy_peaks_from_spectra
# copy_peaks_to_spectra




# %% --------------------- Example 11: peak fitting ----------------------- %% #

fit_peak
fit_peaks

x = np.linspace(0, 10, 100)
y = br.gaussian_fwhm(x=x, amp=4, c=4, w=1) + br.gaussian_fwhm(x=x, amp=2, c=8, w=1)
s = br.Spectrum(x, y)

# fit one peak
s.fit_peak(ranges=(0, 5))
print(s.peaks)
one_peak = s.peaks.calculate_spectrum()

# find all peaks and fit
s.find_peaks()
s.fit_peaks() 
fit = s.peaks.calculate_spectrum()

# quick plot
br.figure()  
s.plot(marker='o') 
fit.plot()
one_peak.plot()   

# %% ----------------- Example 12: defining new methods ------------------- %% #
ss1 = br.Spectra()

# define new function that acts on x and y
def new_function(self, value):
    """Adds a value to the y-axis for each spectra inside Spectra object"""
    for s in self:
        s.y = s.y + value

# attaches function to type Spectrum
br.Spectra.new_function = new_function

# from now one, every spectrum type will have new_function defined.
# even spectrum types created before defying new_function
ss2 = br.Spectrum()

ss1.new_function
ss2.new_function
