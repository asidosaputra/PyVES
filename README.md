## PyVES v1.0.0 : Python library for Vertical Electrical Sounding(VES) modelling

Vertical Electrical Sounding (VES) is by far the most used method for geoelectric surveying, because it is one of the cheapest geophysical method and it gives very good results in many area of interest. The principle of this method is to insert a electric current, of known intensity, through the ground with the help of two electrodes (power electrodes – AB) and measuring the electric potential difference with another two electrodes (measuring electrodes – MN). The investigation depth is proportional with the distance between the power electrodes.

This package support for forward and inverse modelling. Background mathematics of forward modelling by Linear Filtering, this library was provide many of filter coeffcient, such as guptasarma coefficient(7, 11, 22), Das and Koefood filter. Inversion Algortihm are dumped least square inversion, Levenberg-Marquartz(LM) and Single Value Decoompotition(SVD)


## Forward Modelling


## Inversion Modelling
Eath model result of inversion
![inversion](https://github.com/asidosaputra/PyVES/blob/master/examples/Inversion_Earth_Model.png)

Error each iteration
<p align = "center" >
  <img src="https://github.com/asidosaputra/PyVES/blob/master/examples/Inversion_Error.png">
</p>
