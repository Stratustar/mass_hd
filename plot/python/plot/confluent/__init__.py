################################################################################
#
# Plotting routines for the confluent-wet model
#
# Scope is deliberately narrow: these routines assume a confluent-wet frame and
# nothing else. They read its field names directly instead of resolving between
# models, so handing them a frame from another model fails immediately on the
# missing attribute rather than silently plotting the wrong quantity.
#
# Same convention as plot/plot_hd.py and plot/lyotropic.py: every routine takes
# a `frame` and an `engine` (defaults to pyplot) and draws directly onto it.
# Time-series routines take the archive `oa` instead. Nothing here loads
# archives or writes files.
#
################################################################################

from .fields import velocity, velocity_field, velocity_rms

__all__ = ["velocity", "velocity_field", "velocity_rms"]
