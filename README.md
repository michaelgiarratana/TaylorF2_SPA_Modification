February 28, 2025:
------------------

```
bilby.gwpy
```

```
import lal
import gwpy
from gwpy.timeseries import TimeSeries
import bilby

data = TimeSeries.read("/home/giarratana/X-Test-1197008864-16.gwf", channel="Test")
ifo=bilby.gw.detector.get_empty_interferometer("H1")
ifo.strain_data.set_from_gwpy_timeseries(data)
```
