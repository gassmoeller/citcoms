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

from CitcomSLib import CPU_time, output, output_time, output_checkpoint, return_dt, return_t, return_step
from pyre.components.Component import Component
import journal



class Solver(Component):

    def __init__(self, name, facility="solver"):
        Component.__init__(self, name, facility)

        self.all_variables = None
        self.communicator = None
        self.start_cpu_time = 0
        return


    def _dt(self):
        '''get the value of dt from the C code'''
        return return_dt(self.all_variables)


    def _t(self):
        '''get the value of t from the C code'''
        return return_t(self.all_variables)


    def _step(self):
        '''get the value of step from the C code'''
        return return_step(self.all_variables)


    # Set these attributes as read-only properties, so that they are
    # always in accordance with their counterparts in the C code
    t = property(_t)
    dt = property(_dt)
    step = property(_step)


    def initialize(self, application):
        from CitcomSLib import citcom_init, global_default_values, set_signal

        # initialize the "struct All_variables E" in C
        comm = application.solverCommunicator
        all_variables = citcom_init(comm.handle())
        self.communicator = comm
        self.all_variables = all_variables

        # defined in child classes
        self.initializeSolver()

        # information about clock time
        self.start_cpu_time = CPU_time()

        # global interuption handling routine defined once here
        set_signal()

        # default values for various parameters
        global_default_values(self.all_variables)


        # read input parameters from file(s) and command line
        inv = self.inventory

        inv.mesher.initialize(all_variables)
        inv.tsolver.initialize(all_variables)
        inv.vsolver.initialize(all_variables)

        inv.bc.initialize(all_variables)
        inv.const.initialize(all_variables)
        inv.ic.initialize(all_variables)
        inv.output.initialize(all_variables)
        inv.param.initialize(all_variables)
        inv.phase.initialize(all_variables)
        inv.tracer.initialize(all_variables)
        inv.visc.initialize(all_variables)

        from CitcomSLib import return_rank, return_pid
        rank = return_rank(self.all_variables)
        if rank == 0:
            pid = return_pid(self.all_variables)
            stream = open("pid%d.cfg" % pid, "w")
        else:
            stream = None

        # XXX: This is a heck
        # Write controller inventory to pid file
        if stream:
            print >> stream, "[CitcomS.controller]"
            print >> stream, ("monitoringFrequency=%d" %
                application.controller.inventory.monitoringFrequency)
            print >> stream, ("checkpointFrequency=%d" %
                application.controller.inventory.checkpointFrequency)
            print >> stream

        # passing the parameters to C
        self.setProperties(stream)

        self.initial_setup()

        return


    def launch(self, application):
        if self.inventory.ic.inventory.post_p:
            from CitcomSLib import readCheckpoint, postProcessing
            readCheckpoint(self.all_variables)

            # the program will finish after post_processing
            postProcessing(self.all_variables)
            self.save(1)
            self.endSimulation()
            return

        if self.inventory.ic.inventory.restart:
            from CitcomSLib import readCheckpoint
            readCheckpoint(self.all_variables)
        else:
            # initial conditions
            ic = self.inventory.ic
            ic.launch()

            self.solveVelocities()

        # stop the computation if only computes stokes' problem
        if self.inventory.stokes_flow_only:
            self.advectTracers()
            self.save(1)
            self.endSimulation()

        return


    def initial_setup(self):
        mesher = self.inventory.mesher
        vsolver = self.inventory.vsolver
        tsolver = self.inventory.tsolver

        # create mesh
        mesher.run()

        # initialize const. related to mesh
        vsolver.launch()
        tsolver.launch()
        return


    def solveVelocities(self):
        vsolver = self.inventory.vsolver
        vsolver.run()
        return



    def advectTemperature(self, dt):
        tsolver = self.inventory.tsolver
        tsolver.run(dt)
        return



    def advectTracers(self):
        self.inventory.tracer.run()
        return



    def newStep(self):
        return



    def stableTimestep(self):
        tsolver = self.inventory.tsolver
        dt = tsolver.stable_timestep()
        return dt



    def advance(self, dt):
        self.advectTemperature(dt)
        self.advectTracers()
        self.solveVelocities()
        return


    def endTimestep(self, done):
        self.inventory.visc.updateMaterial()
        self.inventory.bc.updatePlateVelocity()
        self.inventory.bc.updatePlateTemperature()
        return done


    def endSimulation(self):
        self._avgCPUTime()
        from CitcomSLib import citcom_finalize
        citcom_finalize(self.all_variables, 0)
        return


    def _avgCPUTime(self):
        step = self.step
        total_cpu_time = CPU_time() - self.start_cpu_time

        rank = self.communicator.rank
        if not rank:
            import sys
            sys.stderr.write("Average cpu time taken for velocity step = %f\n"
                             % (total_cpu_time / (step+1)) )
        return


    def save(self, monitoringFrequency):
        step = self.step

        # output spacing is 'monitoringFrequency'
        if not (step % monitoringFrequency):
            output(self.all_variables, step)

        output_time(self.all_variables, step)
        return


    def checkpoint(self, checkpointFrequency):
        step = self.step

        if not (step % checkpointFrequency):
            output_checkpoint(self.all_variables)
        return


    def setProperties(self, stream):

        from CitcomSLib import Solver_set_properties, check_settings_consistency

        Solver_set_properties(self.all_variables, self.inventory, stream)

        inv = self.inventory
        inv.vsolver.setProperties(stream)
        inv.mesher.setProperties(stream)
        inv.tsolver.setProperties(stream)

        inv.bc.setProperties(stream)
        inv.const.setProperties(stream)
        inv.ic.setProperties(stream)
        inv.output.setProperties(stream)
        inv.param.setProperties(stream)
        inv.phase.setProperties(stream)
        inv.tracer.setProperties(stream)
        inv.visc.setProperties(stream)

        check_settings_consistency(self.all_variables);

        return


    class Inventory(Component.Inventory):

        import pyre.inventory as inv

        # component modules
        import CitcomS.Components.Advection_diffusion as Advection_diffusion
        import CitcomS.Components.Stokes_solver as Stokes_solver

        # components
        from CitcomS.Components.BC import BC
        from CitcomS.Components.Const import Const
        from CitcomS.Components.IC import IC
        from CitcomS.Components.Output import Output
        from CitcomS.Components.Param import Param
        from CitcomS.Components.Phase import Phase
        from CitcomS.Components.Tracer import Tracer
        from CitcomS.Components.Visc import Visc


        tsolver = inv.facility("tsolver",
                               factory=Advection_diffusion.temperature_diffadv)
        vsolver = inv.facility("vsolver",
                               factory=Stokes_solver.incompressibleNewtonian)

        bc = inv.facility("bc", factory=BC)
        const = inv.facility("const", factory=Const)
        ic = inv.facility("ic", factory=IC)
        output = inv.facility("output", factory=Output)
        param = inv.facility("param", factory=Param)
        phase = inv.facility("phase", factory=Phase)
        tracer = inv.facility("tracer", factory=Tracer)
        visc = inv.facility("visc", factory=Visc)

        datadir = inv.str("datadir", default=".")
        datadir_old = inv.str("datadir_old", default=".")

        rayleigh = inv.float("rayleigh", default=1e+05)
        dissipation_number = inv.float("dissipation_number", default=0.0)
        gruneisen = inv.float("gruneisen", default=0.0)
        surfaceT = inv.float("surfaceT", default=0.1)
        adiabaticT0 = inv.float("adiabaticT0", default=0.4)
        Q0 = inv.float("Q0", default=0.0)

        stokes_flow_only = inv.bool("stokes_flow_only", default=False)

        verbose = inv.bool("verbose", default=False)
        see_convergence = inv.bool("see_convergence", default=True)


# version
__id__ = "$Id$"

# End of file
