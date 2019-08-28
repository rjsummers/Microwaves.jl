module Radar

using Unitful
using PhysicalConstants.CODATA2018

c_0 = Unitful.ustrip(SpeedOfLightInVacuum)

function range(return_time::Real)
    """
    Returns the range to a target in meters, given the return time in
    seconds.
    """
    c_0 * return_time / 2
end

function unambiguous_range(prf::Real)
    """
    Returns the maximum unambiguous range in meters given
    a radar system's prf in Hz.
    """
    c_0 / (2 * prf)
end

function range_resolution(pulse_width::Real)
    """
    Returns the range resolution in meters of a radar with a
    pulse width specified in seconds.
    """
    c_0 * pulse_width / 2
end

end # module