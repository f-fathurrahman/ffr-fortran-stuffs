The subroutine should be defined inside a module.

Using emulation mode:

```
pgf90 -Mcuda=emu prog.cuf
```

The attribute `global` indicates that the code is to run on the device but is called
from the host.
The term `global` describes the scope; the subroutine is seen from both the host and
the device.


