# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""H5MD trajectories --- :mod:`MDAnalysis.coordinates.H5MD`
========================================================

The `H5MD`_ trajectory file format is based upon the general, high performance
`HDF5`_ file format.
HDF5 files are self documenting and can be accessed with the `h5py`_ library.
The HDF5 library (and `H5PY`_) must be installed; otherwise, H5MD files
cannot be read by MDAnalysis. If `H5PY`_ is not installed, a ``RuntimeError``
is raised.

HDF5 can make use of parallel file system features through the MPI-IO
interface of the HDF5 library to improve parallel reads and writes.


The `H5MD`_ file format is based upon `HDF5`_, which makes use of parallel
file system features through the MPI-IO interface of the HDF5 library.
The reader currently uses the `H5PY`_ library to access data from an H5MD file.

.. _`H5MD`: https://nongnu.org/h5md/index.html
.. _`HDF5`: https://www.hdfgroup.org/solutions/hdf5/
.. _`H5PY`: http://docs.h5py.org/
.. _`H5MD notation`: https://nongnu.org/h5md/modules/units.html
.. _`MDAnalysis notation`: https://userguide.mdanalysis.org/1.0.0/units.html


Classes
-------

.. autoclass:: Timestep
   :members:

   .. attribute:: positions

      coordinates of the atoms as a :class:`numpy.ndarray` of shape
      `(n_atoms, 3)`

   .. attribute:: velocities

      velocities of the atoms as a :class:`numpy.ndarray` of shape
      `(n_atoms, 3)`; only available if the trajectory contains velocities
      or if the `velocities` = ``True`` keyword has been supplied.

   .. attribute:: forces

      forces of the atoms as a :class:`numpy.ndarray` of shape
      `(n_atoms, 3)`; only available if the trajectory contains forces
      or if the `forces` = ``True`` keyword has been supplied.


.. autoclass:: H5MDReader
   :members:

"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates import base, core
from MDAnalysis.exceptions import NoDataError
try:
    import h5py
except ImportError:
    HAS_H5PY = False
else:
    HAS_H5PY = True


class Timestep(base.Timestep):
    """H5MD Timestep
    """
    order = 'C'

    def _init_unitcell(self):
        return np.zeros((3, 3), dtype=np.float32)

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)

        lengths *A*, *B*, *C* are in the MDAnalysis length unit (Å), and
        angles are in degrees.

        Setting dimensions will populate the underlying native format
        description (triclinic box vectors). If `edges
        <https://nongnu.org/h5md/h5md.html#simulation-box>`_ is a matrix,
        the box is of triclinic shape with the edge vectors given by
        the rows of the matrix.
        """
        return core.triclinic_box(*self._unitcell)

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell[:] = core.triclinic_vectors(box)


class H5MDReader(base.ReaderBase):
    """Reader for the H5MD format.

    See `h5md documentation <https://nongnu.org/h5md/h5md.html>`_
    for a detailed overview of the H5MD file format.

    .. rubric:: Units

    Units are read from the attributes of the position, velocity, force,
    and time datasets. The unit string is translated from `H5MD notation`_ to
    `MDAnalysis notation`_. If MDAnalysis does not recognize the unit
    provided by the H5MD file, a ``RuntimeError`` is raised. If MDAnalysis
    does not recognize the units, it is likely because that unit string is
    not defined in MDAnalysis. If no units are provided,
    MDAnalysis stores a value of ``None`` for each unit. MDAnalysis will raise
    a ``ValueError`` if ``convert_units`` is set to ``True``, yet the
    `H5MD`_ file does not contain units. Set ``convert_units=False`` if
    the `H5MD`_ file contains no units.

    Currently reads .h5md files with the following HDF5 hierarchy:

    .. code-block:: text

        Notation:
        (name) is an HDF5 group that the reader recognizes
        {name} is an HDF5 group with arbitrary name
        [variable] is an HDF5 dataset
        <dtype> is dataset datatype
        +-- is an attribute of a group or dataset

        H5MD root
         \-- (h5md)
            +-- version <int>
            \-- author
                +-- name <str>, author's name
                +-- email <str>, optional email address
            \-- creator
                +-- name <str>, file that created .h5md file
                +-- version
         \-- (particles)
            \-- {group1}
                \-- (box)
                    +-- dimension : <int>, number of spatial dimensions
                    +-- boundary : <str>, boundary conditions of unit cell
                    \-- (edges)
                        \-- [step] <int>, gives frame
                        \-- [value] <float>, gives box dimensions
                \-- (position)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- units <str>
                    \-- [value] <float>, gives numpy arrary of positions
                                         with shape (n_atoms, 3)
                        +-- units <str>
                \-- (velocity)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- units <str>
                    \-- [value] <float>, gives numpy arrary of velocities
                                         with shape (n_atoms, 3)
                        +-- units <str>
                \-- (force)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- units <str>
                    \-- [value] <float>, gives numpy arrary of forces
                                         with shape (n_atoms, 3)
                        +-- units <str>
         \-- (observables)
            \-- (lambda)
                \-- [step] <int>, gives frame
                \-- [time] <float>, gives time
                \-- [value] <float>
            \-- (step)
                \-- [step] <int>, gives frame
                \-- [time] <float>, gives time
                \-- [value] <int>, gives integration step

    .. note::
        The reader does not currently read mass or charge data.

    .. versionadded:: 2.0.0

    """

    format = 'H5MD'
    # units is defined as instance-level variable and set from the
    # H5MD file in __init__() below

    # This dictionary is used to translate H5MD units to MDAnalysis units.
    # (https://nongnu.org/h5md/modules/units.html)
    _unit_translation = {
        'time': {
            'ps': 'ps',
            'fs': 'fs',
            'ns': 'ns',
            'second': 's',
            'sec': 's',
            's': 's',
            'AKMA': 'AKMA',
            },
        'length': {
            'Angstrom': 'Angstrom',
            'angstrom': 'Angstrom',
            'A': 'Angstrom',
            'nm': 'nm',
            'pm': 'pm',
            'fm': 'fm',
            },
        'velocity': {
            'Angstrom ps-1': 'Angstrom/ps',
            'A ps-1': 'Angstrom/ps',
            'Angstrom fs-1': 'Angstrom/fs',
            'A fs-1': 'Angstrom/fs',
            'Angstrom AKMA-1': 'Angstrom/AKMA',
            'A AKMA-1': 'Angstrom/AKMA',
            'nm ps-1': 'nm/ps',
            'nm ns-1': 'nm/ns',
            'pm ps-1': 'pm/ps',
            'm s-1': 'm/s'
            },
        'force':  {
            'kJ mol-1 Angstrom-1': 'kJ/(mol*Angstrom)',
            'kJ mol-1 nm-1': 'kJ/(mol*nm)',
            'Newton': 'Newton',
            'N': 'N',
            'J m-1': 'J/m',
            'kcal mol-1 Angstrom-1': 'kcal/(mol*Angstrom)',
            'kcal mol-1 A-1': 'kcal/(mol*Angstrom)'
            }
    }
    _Timestep = Timestep
    _data_keywords = ('step', 'lambda')

    def __init__(self, filename, convert_units=True, **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`h5py.File`
            trajectory filename or open h5py file
        convert_units : bool (optional)
            convert units to MDAnalysis units
        **kwargs : dict
            General reader arguments.

        Raises
        ------
        RuntimeError
            when `H5PY`_ is not installed
        RuntimeError
            when a unit is not recognized by MDAnalysis
        ValueError
            when ``n_atoms`` changes values between timesteps
        ValueError
            when ``convert_units=True`` but the H5MD file contains no units
        NoDataError
            when the H5MD file has no 'position', 'velocity', or
            'force' group

        """
        if not HAS_H5PY:
            raise RuntimeError("Please install h5py")
        super(H5MDReader, self).__init__(filename, **kwargs)
        self.filename = filename
        ## add parallel MPI kwwargs
        ## kwwargs needs to get pulled out
        self.open_trajectory()

        # _has dictionary used for checking whether h5md file has
        # 'position', 'velocity', or 'force' groups in the file
        self._has = {name: name in self._particle_group for
                     name in ('position', 'velocity', 'force')}

        # Gets n_atoms from first available group
        for name, value in self._has.items():
            if value:
                self.n_atoms = self._particle_group[name]['value'].shape[1]
                break
        else:
            raise NoDataError("Provide at least a position, velocity"
                              " or force group in the h5md file.")
        self.ts = self._Timestep(self.n_atoms,
                                 velocities=self.has_velocities,
                                 forces=self.has_forces,
                                 **self._ts_kwargs)
        self.units = {'time': None,
                      'length': None,
                      'velocity': None,
                      'force': None}
        self._set_translated_units()  # fills units dictionary
        self._read_next_timestep()

    @staticmethod
    def _format_hint(thing):
        """Can this Reader read *thing*"""
        # nb, filename strings can still get passed through if
        # format='H5MD' is used
        return HAS_H5PY and isinstance(thing, h5py.File)

    def open_trajectory(self):
        """opens the trajectory file using h5py library"""
        self._frame = -1
        if isinstance(self.filename, h5py.File):
            self._file = self.filename
        else:
            self._file = h5py.File(self.filename, 'r')
        # pulls first key out of 'particles'
        # allows for arbitrary name of group1 in 'particles'
        self._particle_group = self._file['particles'][
                               list(self._file['particles'])[0]]

    def close(self):
        """close reader"""
        self._file.close()
        ## need to change this so it doesnt close the file, since
        ## file might not be reopened in the same away

    @property
    def has_positions(self):
        """True if 'position' group is in trajectory."""
        return self._has['position']

    @has_positions.setter
    def has_positions(self, value: bool):
        self._has['position'] = value

    @property
    def has_velocities(self):
        """True if 'velocity' group is in trajectory."""
        return self._has['velocity']

    @has_velocities.setter
    def has_velocities(self, value: bool):
        self._has['velocity'] = value

    @property
    def has_forces(self):
        """True if 'force' group is in trajectory."""
        return self._has['force']

    @has_forces.setter
    def has_forces(self, value: bool):
        self._has['force'] = value

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        for name, value in self._has.items():
            if value:
                return self._particle_group[name]['value'].shape[0]

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):
        """reads data from h5md file and copies to current timestep"""
        try:
            for name, value in self._has.items():
                if value:
                    myframe = self._particle_group[name]['step'][frame]
                    break
            else:
                raise NoDataError("Provide at least a position, velocity"
                                  " or force group in the h5md file.")
        except ValueError:
            raise IOError from None

        self._frame = frame
        ts = self.ts
        particle_group = self._particle_group
        ts.frame = frame

        # this block populates the data dictionary
        # DT not read
        if 'observables' in self._file:
            data = self._file['observables']
            for name in self._data_keywords:
                self._copy_data(name, data)
        for name, value in self._has.items():
            if value:
                self.ts.data['time'] = particle_group[name]['time'][frame]
                break

        # set frame box dimensions
        # set triclinic box vectors
        ts._unitcell[:] = particle_group['box/edges/value'][frame, :]

        # set the timestep positions, velocities, and forces with
        # set_ts_attribute() method
        ts.positions = self._set_ts_attribute('position')
        ts.velocities = self._set_ts_attribute('velocity')
        ts.forces = self._set_ts_attribute('force')

        # unit conversion
        if self.convert_units:
            # ensures h5md file has units if convert_units=True
            self._check_and_convert_units()

        return ts

    def _read_next_timestep(self):
        """read next frame in trajectory"""
        return self._read_frame(self._frame + 1)

    def _set_ts_attribute(self, dataset):
        if self._has[dataset]:
            frame_dataset = self._particle_group[
                            dataset]['value'][self.frame, :]
            n_atoms_now = frame_dataset.shape[0]
            if n_atoms_now != self.n_atoms:
                raise ValueError("Frame {} has {} atoms but the initial frame"
                                 " has {} atoms. MDAnalysis is unable to deal"
                                 " with variable topology!"
                                 "".format(self.frame,
                                           n_atoms_now,
                                           self.n_atoms))
            else:
                return frame_dataset

    def _check_and_convert_units(self):

        errmsg = "H5MD file must have readable units if ``convert_units`` is"
        " set to ``True`` MDAnalysis sets ``convert_units=True`` by default."
        " Set ``convert_units=False`` to to load Universe."

        if self.units['time'] == None:
            raise ValueError(errmsg)
        else:
            self.ts.time = self.convert_time_from_native(self.ts.time)

        if self._has['position']:
            if self.units['length'] == None:
                raise ValueError(errmsg)
            else:
                self.convert_pos_from_native(self.ts.positions)
                self.convert_pos_from_native(self.ts.dimensions[:3])

        if self._has['velocity']:
            if self.units['velocity'] == None:
                raise ValueError(errmsg)
            else:
                self.convert_velocities_from_native(self.ts.velocities)

        if self._has['force']:
            if self.units['force'] == None:
                raise ValueError(errmsg)
            else:
                self.convert_forces_from_native(self.ts.forces)


    def _set_translated_units(self):
        """converts units from H5MD to MDAnalysis notation
        and fills units dictionary"""

        # need this dictionary to associate 'position': 'length'
        _group_unit_dict = {'time': 'time',
                            'position': 'length',
                            'velocity': 'velocity',
                            'force': 'force'
                            }

        for group, unit in _group_unit_dict.items():
            self._translate_h5md_units(group, unit)

    def _translate_h5md_units(self, group, unit):
        """stores the translated unit string into the units dictionary"""

        errmsg = "{} unit '{}' is not recognized by H5MDReader. Please raise"
        " an issue in https://github.com/MDAnalysis/mdanalysis/issues"

        # doing time unit separately because time has to fish for
        # first available parent group, either position, velocity, or force
        if unit == 'time':
            for name, value in self._has.items():
                if value:
                    if 'units' in self._particle_group[name]['time'].attrs:
                        try:
                            self.units['time'] = self._unit_translation[
                            'time'][self._particle_group[name][
                            'time'].attrs['units']]
                            break
                        except KeyError:
                            raise RuntimeError(errmsg.format(unit,
                                               self._particle_group[
                                               name]['time'].attrs[
                                               'units'])) from None

        else:
            if self._has[group]:
                if 'units' in self._particle_group[group].attrs:
                    try:
                        self.units[unit] = self._unit_translation[unit][
                        self._particle_group[group].attrs['units']]
                    except KeyError:
                        raise RuntimeError(errmsg.format(unit,
                                           self._particle_group[group].attrs[
                                           'units'])) from None

    def _copy_data(self, name, group):
        """assigns values to keys in data dictionary"""
        self.ts.data[name] = group[name]['value'][self.frame]


class H5MDWriter(base.WriterBase):
    """Writer for the H5MD format.

    Classes
    -------

    .. autoclass:: Timestep
       :members:

       .. attribute:: positions

          coordinates of the atoms as a :class:`numpy.ndarray` of shape
          `(n_atoms, 3)`; only available if the trajectory contains positions
          or if the ``positions = True`` keyword has been supplied.

       .. attribute:: velocities

          velocities of the atoms as a :class:`numpy.ndarray` of shape
          `(n_atoms, 3)`; only available if the trajectory contains velocities
          or if the ``velocities = True`` keyword has been supplied.

       .. attribute:: forces

          forces of the atoms as a :class:`numpy.ndarray` of shape
          `(n_atoms, 3)`; only available if the trajectory contains forces
          or if the ``forces = True`` keyword has been supplied.


    .. autoclass:: H5MDReader
       :members:

    .. versionadded:: 2.0.0

    """

    format = 'H5MD'
    multiframe = True
    _unit_translation = {
        'time': {
            'ps': 'ps',
            'fs': 'fs',
            'ns': 'ns',
            'second': 's',
            'sec': 's',
            's': 's',
            'AKMA': 'AKMA',
            },
        'length': {
            'Angstrom': 'Angstrom',
            'angstrom': 'Angstrom',
            'A': 'Angstrom',
            'nm': 'nm',
            'pm': 'pm',
            'fm': 'fm',
            },
        'velocity': {
            'Angstrom/ps': 'Angstrom ps-1',
            'A/ps': 'Angstrom ps-1',
            'Angstrom/fs': 'Angstrom fs-1',
            'A/fs': 'Angstrom fs-1',
            'Angstrom/AKMA': 'Angstrom AKMA-1',
            'A/AKMA': 'Angstrom AKMA-1',
            'nm/ps': 'nm ps-1',
            'nm/ns': 'nm ns-1',
            'pm/ps': 'pm ps-1',
            'm/s': 'm s-1'
            },
        'force':  {
            'kJ/(mol*Angstrom)': 'kJ mol-1 Angstrom-1',
            'kJ/(mol*nm)': 'kJ mol-1 nm-1',
            'Newton': 'Newton',
            'N': 'N',
            'J/m': 'J m-1',
            'kcal/(mol*Angstrom)': 'kcal mol-1 Angstrom-1',
            'kcal/(mol*A)': 'kcal mol-1 Angstrom-1'
            }
    }

    def __init__(self,
                 filename,
                 n_atoms,
                 convert_units=True,
                 units=None,
                 **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`h5py.File`
            trajectory filename or open h5py file
        n_atoms : int
            number of atoms in trajectory
        convert_units : bool (optional)
            convert units to MDAnalysis units
        **kwargs : dict

        """
        self.filename = filename
        if n_atoms == 0:
            raise ValueError("H5MDWriter: no atoms in output trajectory")
        self.n_atoms = n_atoms
        self.convert_units = convert_units

        self.h5mdfile = None
        self.has_positions = kwargs.get('positions', False)
        self.has_velocities = kwargs.get('velocities', False)
        self.has_forces = kwargs.get('forces', False)

        # Pull out various keywords to store in 'h5md' group
        self.author = kwargs.pop('author', 'N/A')
        self.author_email = kwargs.pop('author_email', None)
        self.creator = kwargs.pop('creator', 'MDAnalysis')
        self.creator_version = kwargs.pop('creator_version', mda.__version__)

        # Units stored as instance variable
        # Will be filled after translated from ts object
        self.units = units

    def _init_h5md(self, periodic=True):
        """Initializes H5MD trajectory.

        The `H5MD`_ file is opened using the `H5PY`_ library. """

        h5md_file = h5py.File(self.filename, 'w')
        h5md_group = h5md_file.create_group('h5md')
        h5md_group.attrs['version'] = np.array([1,1])
        author_group = h5md_file.create_group('h5md/author')
        author_group.attrs['name'] = self.author
        if self.author_email is not None:
            author_group.attrs['email'] = self.author_email
        creator_group = h5md_file.create_group('h5md/creator')
        creator_group.attrs['name'] = self.creator
        creator_group.attrs['version'] = self.creator_version

        self.h5mdfile = h5md_file
        self._trajectory_group = self._create_particles_group(name='trajectory')

    def is_periodic(self, ts):
        """Test if timestep ``ts`` contains a periodic box.

        Parameters
        ----------
        ts : :class:`Timestep`
             :class:`Timestep` instance containing coordinates to
             be written to trajectory file

        Returns
        -------
        bool
            Return ``True`` if `ts` contains a valid simulation box
        """
        return np.all(ts.dimensions > 0)

    def _write_next_frame(self, ag):
        """Write information associated with ``ag`` at current frame
        into trajectory

        Parameters
        ----------
        ag : AtomGroup or Universe

        .. versionadded:: 2.0.0

        """
        try:
            # Atomgroup?
            ts = ag.ts
        except AttributeError:
            try:
                # Universe?
                ts = ag.trajectory.ts
            except AttributeError:
                errmsg = "Input obj is neither an AtomGroup or Universe"
                raise TypeError(errmsg) from None

        if ts.n_atoms != self.n_atoms:
            raise IOError("H5MDWriter: Timestep does not have"
                          " the correct number of atoms")

        # Opens H5MD file to write with H5PY library and initializes
        # 'h5md' group and checks if unitcell is periodic
        if self.h5mdfile is None:
            self._init_h5md(periodic = self.is_periodic(ts))

        trajectory = self._trajectory_group  # conventient assignment

        # Checks to see if writer given boolean arguments for position,
        # velocity, and force, then instantiates H5MD dataset with
        # _create_element() method. Also adds translated units to each group.
        # Box dimensions are created at first available boolean
        if self.has_positions:
            self.position = self._create_element(parent_group=trajectory,
                                                 name='position',
                                                 data=ts.positions,
                                                 dataset_type='time')
            trajectory['position'].attrs['units'] = self._unit_translation['length'][self.units['length']]
            self.box = self._create_box(dimension=3,
                                        boundary=['periodic',
                                                  'periodic',
                                                  'periodic'],
                                        data=ts.triclinic_dimensions,
                                        dataset_type='time',
                                        step_from=self.position)

        if self.has_velocities:
            if self.has_positions:
                self.velocity = self._create_element(parent_group=trajectory,
                                                     name='velocity',
                                                     data=ts.velocities,
                                                     step_from=self.position,
                                                     dataset_type='time')
                trajectory['velocity'].attrs['units'] = self._unit_translation['velocity'][self.units['velocity']]
            else:
                self.velocity = self._create_element(parent_group=trajectory,
                                                     name='velocity',
                                                     data=ts.velocities,
                                                     dataset_type='time')
                trajectory['velocity'].attrs['units'] = self._unit_translation['velocity'][self.units['velocity']]
                self.box = self._create_box(dimension=3,
                                            boundary=['periodic',
                                                      'periodic',
                                                      'periodic'],
                                            data=ts.triclinic_dimensions,
                                            step_from=self.velocity)

        if self.has_forces:
            if self.has_positions:
                self.force = self._create_element(parent_group=trajectory,
                                                  name='force',
                                                  data=ts.forces,
                                                  step_from=self.position,
                                                  dataset_type='time')
                trajectory['force'].attrs['units'] = self._unit_translation['force'][self.units['force']]
            elif self.has_velocities:
                self.force = self._create_element(parent_group=trajectory,
                                                  name='force',
                                                  data=ts.forces,
                                                  step_from=self.velocity,
                                                  dataset_type='time')
                trajectory['force'].attrs['units'] = self._unit_translation['force'][self.units['force']]
            else:
                self.force = self._create_element(parent_group=trajectory,
                                                  name='force',
                                                  data=ts.forces,
                                                  dataset_type='time')
                trajectory['force'].attrs['units'] = self._unit_translation['force'][self.units['force']]
                self.box = self._create_box(dimension=3,
                                            boundary=['periodic',
                                                      'periodic',
                                                      'periodic'],
                                            data=ts.triclinic_dimensions,
                                            step_from=self.force)


        return self._write_next_timestep(ts)

    def _write_next_timestep(self, ts):
        """Write coordinates and unitcell information to H5MD file.

        Do not call this method directly; instead use
        :meth:`write` because some essential setup is done
        there before writing the first frame.

        """
        pos = ts.positions
        vel = ts.velocities
        force = ts.forces
        time = ts.time
        unitcell = ts.dimensions

        #if self.convert_units:
            #pos = self.convert_pos_to_native(pos, inplace=False)
            #vel = self.convert_velocities_to_native(vel, inplace=False)
            #force = self.convert_forces_to_native(force, inplace=False)
            #time = self.convert_time_to_native(time, inplace=False)
            #unitcell = self.convert_dimensions_to_unitcell(ts)

        if self.has_positions:
            self.position.append(value=ts.positions,
                                 step=ts.frame,
                                 time=ts.time)
        if self.has_velocities:
            self.velocity.append(value=ts.velocities,
                                 step=ts.frame,
                                 time=ts.time)
        if self.has_forces:
            self.force.append(value=ts.forces,
                              step=ts.frame,
                              time=ts.time)
        self.box.edges.append(value=ts.triclinic_dimensions,
                              step=ts.frame,
                              time=ts.time)


    def _create_particles_group(self, name):
        """ """
        return self.h5mdfile.require_group('particles').require_group(name)

############################ PIERRE CODE BELOW ##############################################
    def _create_box(self, dimension=None, boundary=None, **kwargs):
        """ """
        self.box = self._trajectory_group.require_group('box')
        if dimension is not None and boundary is not None:
            assert len(boundary)==dimension
            self.box.attrs['dimension'] = dimension
            self.box.attrs['boundary'] = np.string_(boundary)
        if len(kwargs)>0:

            self.box.edges = self._create_element(parent_group=self.box, name='edges', **kwargs)

        return self.box

    def _create_element(self, parent_group, name, **kwargs):
        """ """
        if name in parent_group:
            tmp_element = parent_group[name]
            assert 'value' in tmp_element
            assert 'step' in tmp_element

            return Element(parent_group, name, **kwargs)

        data = kwargs.pop('data', None)
        if data is not None:
            if type(data) is not np.ndarray:
                data = np.asarray(data)
            kwargs['shape'] = (0,) + data.shape
            if 'maxshape' in kwargs:
                kwargs['maxshape'] = (None,) + kwargs['maxshape']
            else:
                kwargs['maxshape'] = (None,) + data.shape
            kwargs['dtype'] = data.dtype
        return Element(parent_group, name, **kwargs)

    def _default_chunks(shape):
        """ """
        CHUNK_0_MAX=32
        CHUNK_1_MAX=128
        result = list(shape)
        if len(shape)==0:
            raise ValueError
        if shape==(0,):
            return (64,)
        if shape[0]==0:
            result[0] = CHUNK_0_MAX
            result[1] = min(CHUNK_1_MAX, shape[1])
            return tuple(result)
        else:
            result[:2] = CHUNK_0_MAX, CHUNK_1_MAX
            return tuple(np.minimum(result,shape))

    def close(self):
        """ """
        if self.h5mdfile is not None:
            self.h5mdfile.close()
            self.h5mdfile = None

class Element(h5py.Group):
    """H5MD element"""


    def __init__(self, parent_group, name, **kwargs):
        """ Documenation
        """

        self.dataset_type = kwargs.pop('dataset_type', None)
        self.own_step=False

        if self.dataset_type == 'time':
            is_new = name not in parent_group
            g = parent_group.require_group(name)
            if is_new:
                step_from = kwargs.pop('step_from', None)
                if step_from is not None:
                    g['step'] = step_from.step
                    self.step = step_from.step
                    self.own_step = False
                else:
                    self.step = g.create_dataset('step', dtype=int, shape=(0,), maxshape=(None,))
                    self.own_step = True
                if self.own_step:
                    self.time = g.create_dataset('time', dtype=float, shape=(0,), maxshape=(None,))
                else:
                    g['time'] = self.time = step_from.time
                self.value = g.create_dataset('value', **kwargs)
            else:
                self.step = g['step']
                if 'offset' in g['step'].attrs:
                    self.step_offset = g['step'].attrs['offset']
                else:
                    self.step_offset = None
                self.time = g['time']
                if 'offset' in g['time'].attrs:
                    self.time_offset = g['time'].attrs['offset']
                else:
                    self.time_offset = None
                self.value = g['value']
                super(Element, self).__init__(g._id)

        elif self.dataset_type == 'fixed':
            if name in parent_group:
                dset = parent_group[name]
            else:
                dset = parent_group.create_dataset(name, **kwargs)
            super(Element, self).__init__(dset._id)
            self.step = None
            self.step_offset = None
            self.time = None
            self.time_offset = None

    def append(self, value, step, time=None, region=None, collective=False):
        """ """
        if self.dataset_type == 'time':
            self.value.resize(self.value.shape[0]+1, axis=0)
            if region is not None:
                if collective:
                    with self.value.collective:
                        self.value[-1,region[0]:region[1],...] = value
                else:
                    self.value[-1,region[0]:region[1],...] = value
            else:
                self.value[-1] = value
            if self.own_step:
                self.step.resize(self.step.shape[0]+1, axis=0)
                self.step[-1] = step
                if self.time and len(self.time.shape)==1:
                    self.time.resize(self.time.shape[0]+1, axis=0)
                    self.time[-1] = time

        elif self.dataset_type == 'fixed':
            pass

    def get_by_idx(self, idx):
        return self['value'][idx]

    @property
    def element_type(self):
        if dataset_type == 'time':
            return 'TimeElement'
        elif dataset_type == 'fixed':
            return 'FixedElement'

    def __repr__(self):
        if dataset_type == 'time':
            return 'H5MD TimeElement'
        elif dataset_type == 'fixed':
            return 'H5MD FixedElement'
