#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#<LicenseText>
#
# CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
# Copyright (C) 2002-2005, California Institute of Technology.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#</LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from CitcomComponent import CitcomComponent


class Tracer(CitcomComponent):


    def __init__(self, name="tracer", facility="tracer"):
        CitcomComponent.__init__(self, name, facility)
        return



    def run(self):
        from CitcomSLib import Tracer_tracer_advection
        Tracer_tracer_advection(self.all_variables)
        return



    def setProperties(self, stream):
        # convert lists of strings to lists of ints/floats
        inv = self.inventory
        inv.z_interface = map(float, inv.z_interface)
        inv.buoyancy_ratio = map(float, inv.buoyancy_ratio)
        inv.initial_content = map(float, inv.initial_content)
        inv.Q0_enriched = map(float, inv.Q0_enriched)

        from CitcomSLib import Tracer_set_properties
        Tracer_set_properties(self.all_variables, self.inventory, stream)
        return



    class Inventory(CitcomComponent.Inventory):


        import pyre.inventory as inv


        tracer = inv.bool("tracer", default=False)

        # tracer_ic_method=0 (random generated array)
        # tracer_ic_method=1 (all proc read the same file)
        # tracer_ic_method=2 (each proc reads its own file)
        tracer_ic_method = inv.int("tracer_ic_method", default=0)

        # (tracer_ic_method == 0)
        tracers_per_element = inv.int("tracers_per_element", default=10)

        # (tracer_ic_method == 1)
        tracer_file = inv.str("tracer_file", default="tracer.dat")

        # How many flavors of tracers
        # If tracer_flavors > 0, each element will report the number of
        # tracers of each flavor inside it. This information can be used
        # later for many purposes. One of it is to compute composition,
        # either using absolute method or ratio method.
        tracer_flavors = inv.int("tracer_flavors", default=0)

        # How to initialize tracer flavors
        ic_method_for_flavors = inv.int("ic_method_for_flavors", default=0)
        z_interface = inv.list("z_interface", default=[0.7])
        ictracer_grd_file = inv.str("ictracer_grd_file", default="")
        ictracer_grd_layers = inv.int("ictracer_grd_layers", default=2)

        # Warning level
        itracer_warnings = inv.bool("itracer_warnings", default=True)

        # Enriched internal heat production
        tracer_enriched = inv.bool("tracer_enriched", default=False)
        Q0_enriched = inv.list("Q0_enriched", default=[0.0])

        # Regular grid parameters
        regular_grid_deltheta = inv.float("regular_grid_deltheta", default=1.0)
        regular_grid_delphi = inv.float("regular_grid_delphi", default=1.0)

        # Analytical Test Function
        #analytical_tracer_test = inv.int("analytical_tracer_test", default=0)

        chemical_buoyancy = inv.bool("chemical_buoyancy", default=True)

        # ibuoy_type=0 (absolute method, not implemented)
        # ibuoy_type=1 (ratio method)
        buoy_type = inv.int("buoy_type", default=1)
        buoyancy_ratio = inv.list("buoyancy_ratio", default=[1.0])
        initial_content = inv.list("initial_content", default=[1.0])
        zdep_buoyancy = inv.int("zdep_buoyancy", default=0)
        pressure_oversampling = inv.int("pressure_oversampling", default=1)
        density_file = inv.str("density_file", default="density.dat")
        tdep_buoyancy = inv.int("tdep_buoyancy", default=0)
        start_temp = inv.float("start_temp", default=0)
        end_temp = inv.float("end_temp", default=0)
        ntdeps = inv.int("ntdeps", default=100)

        continents = inv.int("continents", default=0)
        tracer_origin = inv.int("tracer_origin", default=0)
        tracer_origin_set_time = inv.float("tracer_origin_set_time", default=0.0)
        chemical_changes = inv.int("chemical_changes", default=0)
        hotspot_tracks = inv.int("hotspot_tracks", default=0)
        
        # This is not used anymore and is left here for backward compatibility
        reset_initial_composition = inv.bool("reset_initial_composition",
                                             default=False)


        # compositional_rheology=1 (not implemented in this version, TODO:remove)
        #compositional_rheology = inv.bool("compositional_rheology",
        #                                  default=False)
        #compositional_prefactor = inv.float("compositional_prefactor",
        #                                    default=1.0)



# version
__id__ = "$Id$"

# End of file
