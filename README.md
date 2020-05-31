## PyVES v0.0.1 : Python library for Vertical Electrical Sounding(VES) modelling

Vertical Electrical Sounding (VES) is by far the most used method for geoelectric surveying, because it is one of the cheapest geophysical method and it gives very good results in many area of interest. The principle of this method is to insert a electric current, of known intensity, through the ground with the help of two electrodes (power electrodes – AB) and measuring the electric potential difference with another two electrodes (measuring electrodes – MN). The investigation depth is proportional with the distance between the power electrodes.

This package support for forward and inverse modelling. Background mathematics of forward modelling by Linear Filtering(....), this library was provide many of filter coeffcient, such as guptasarma coefficient(7, 11, 22), Das and Koefood filter. Inversion Algortihm are dumped least square inversion, Levenberg-Marquartz(LM) and Single Value Decoompotition(SVD)


## Forward Modelling
<img src="https://github.com/asidosaputra/PyVES/blob/master/examples/Forward_Model.png">

## Inversion Modelling
