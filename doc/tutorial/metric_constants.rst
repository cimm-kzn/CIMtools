Metric constants conversion
===========================

CIMtools\.metric_constants subpackage includes values of masses, volumes, etc for converting to SI values.

Usage is very simple::

    from CIMtools.metric_constants import *

    # convert pounds into kilograms
    2 * lb

    # reverse convert SI into pounds
    2 / lb

    # temperature conversion
    25 * C  # now T in Kelvin
    298 / F  # Kelvin (SI) into Fahrenheit

According to order of math operation, power has higher priority. For prevent errors use ()::

    2 * C ** 2  # wrong
    (2 * C)  ** 2  # first convert, after square
    2 * C * 3  # right


(C) Dr. Ramil Nugmanov
