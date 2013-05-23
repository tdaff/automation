#!/usr/bin/env python

"""
Sqlite backend to save generated structures in a database.

"""

import hashlib
from logging import info, error, debug

from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, Float, String, Text, ForeignKey
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.ext.declarative import declarative_base


# Status codes for if structures are being calculated
NEW = 0
STARTED = 1
FINISHED = 2

Base = declarative_base()

# A secondary table is needed to link the structures to functionalisation sets
sym_associations = Table('sym_associations',
                         Base.metadata,
                         Column('sym_functionalised_structure_id', Integer,
                                ForeignKey('sym_functionalised_structures.id')),
                         Column('functionalisation_id', Integer,
                                ForeignKey('functionalisations.id')))

# keep track of which freeform structures have which groups
freeform_groups = Table('free_associations',
                         Base.metadata,
                         Column('free_functionalised_structure_id', Integer,
                                ForeignKey('free_functionalised_structures.id')),
                         Column('functional_group_id', String,
                                ForeignKey('functional_groups.id')))

# Define models as in the SQLAlchemy tutorial

class SymFunctionalisedStructure(Base):
    """
    SQLAlchemy model for structures that have been derived with the symmetry
    of the base structure. Each .base_structure .name (or .fullname) should
    be the same structure.

    """
    __tablename__ = 'sym_functionalised_structures'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    base_structure = Column(String)
    cif_file = Column(Text)
    status = Column(Integer)

    # Many to many with bidirectional realtionship
    functionalisations = relationship('Functionalisation',
                                      secondary=sym_associations,
                                      backref='sym_functionalised_structures')
    # with the functionalisations relationships we can recreate name:
    # ".".join(["%s@%s" % (func.functional_group_id, func.site)
    #           for func in self.functionalisations]))

    def __init__(self, base_structure, name, cif_file):
        self.base_structure = base_structure
        self.name = name
        self.cif_file = cif_file
        self.status = NEW

    def __repr__(self):
        return "<SymFunctionalisedStructure('[%s]')>" % (self.fullname)

    @property
    def fullname(self):
        """Return the unique composite name"""
        return "%s_func_%s" % (self.base_structure, self.name)


# Functional groups are enumerated separately from the name as the count will
# be variable by mof and functionalisation.

class FunctionalGroup(Base):
    """Essentially the names and identifiers of available groups"""
    __tablename__ = 'functional_groups'

    # Each identifier should be a unique string.
    id = Column(String, primary_key=True, unique=True)
    name = Column(String)

    # The functionalisations should be able to trace back to these
    functionalisations = relationship("Functionalisation",
                                      backref='functional_groups')

    def __init__(self, id, name):
        self.id = id
        self.name = name

    def __repr__(self):
        return "<FunctionalGroup('%s','%s')>" % (self.id, self.name)


class Functionalisation(Base):
    """
    Pairs of functional groups and attachment sites. Note that sites
    with the same name will have no relationship between structures.

    """
    __tablename__ = 'functionalisations'

    id = Column(Integer, primary_key=True)
    site = Column(String)
    functional_group_id = Column(String, ForeignKey('functional_groups.id'))

    def __init__(self, functional_group_id, site):
        self.functional_group_id = functional_group_id
        self.site = site

    def __repr__(self):
        return ("<Functionalisation('%s','%s')>" %
                (self.functional_group_id, self.site))

    def __str__(self):
        return "%s@%s" % (self.functional_group_id, self.site)


# These are not constrained by the symmetry; can be defined by just a string

class FreeFunctionalisedStructure(Base):
    """
    SQLAlchemy model for structures that have freeform functionalisation
    over all the sites of the base structure. Each .base_structure .name
    (or .fullname) should be the same structure but rotations may change
    each generation.

    """
    __tablename__ = 'free_functionalised_structures'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    base_structure = Column(String)
    cif_file = Column(Text)
    status = Column(Integer)

    functional_groups = relationship('FunctionalGroup',
                                     secondary=freeform_groups,
                                     backref='free_functionalised_structures')

    def __init__(self, base_structure, name, cif_file):
        self.base_structure = base_structure
        self.name = name.strip("{}")
        self.cif_file = cif_file
        self.status = NEW

    def __repr__(self):
        return "<FreeFunctionalisedStructure('{%s}')>" % (self.name)

    @property
    def unique_name(self):
        """Return the unique composite name"""
        func_repr = self.name.split(".")
        unique_name = hashlib.md5(str(func_repr)).hexdigest()
        return "%s_free_%s" % (self.base_structure, unique_name)



# Objects to interact with the database in a pluggable way; methods
# can hide the add or create and cacheing

class AlchemyBackend(object):
    """Abstraction for fapswitch+sql interface."""

    def __init__(self, db_name):
        # TODO(tdaff): database is just sqlite the base name for now
        # Eventually this will all be moresql
        self.engine = create_engine('sqlite:///%s.db' % db_name)  #, echo=True)
        # Always create the tables as we often start from a bare structure
        Base.metadata.create_all(self.engine)
        # Bind a session instance to the instance -- easier to query itself
        Session = sessionmaker(bind=self.engine)  # Might need autoflush/commit
        self.session = Session()

    def populate_groups(self, groups):
        """Initialise or update the available groups table from the library."""
        session_groups = [group.id for group in
                          self.session.query(FunctionalGroup)]
        # Assume we just get passed the dict
        for group in groups:
            if not group in session_groups:
                self.session.add(FunctionalGroup(group, groups[group].name))
        self.session.commit()
        self.groups = self.session.query(FunctionalGroup)

    def add_symmetry_structure(self, base_structure, functions, cif_file):
        """
        Insert a structure with symmetry based functionalisation into the
        database. Will try not to duplicate the "group@site" association
        designations.

        """

        existing_func = self.session.query(Functionalisation)

        new_mof_name = ".".join(["@".join(x) for x in functions])
        cif_file = "".join(cif_file)

        # Get an exisiting vesion of the structure if it exists
        structure = self.session.query(SymFunctionalisedStructure).\
            filter(SymFunctionalisedStructure.base_structure == base_structure).\
            filter(SymFunctionalisedStructure.name == new_mof_name).\
            first()

        if structure is None:
            # Not in the database so initialise the structure and then add in
            # the functionalisations
            structure = SymFunctionalisedStructure(base_structure,
                                                   new_mof_name,
                                                   cif_file)

            for group, site in functions:
                # Query for one if it is already stored
                functionalisation = existing_func.\
                    filter(Functionalisation.site == site).\
                    filter(Functionalisation.functional_group_id == group).\
                    first()
                if functionalisation is None:
                    # Create it
                    functionalisation = Functionalisation(group, site)
                    self.session.add(functionalisation)

                # The relationship keeps track of the list
                structure.functionalisations.append(functionalisation)

            self.session.add(structure)
            # New structure inserted
        else:
            structure.cif_file = cif_file

        self.session.commit()

    def add_freeform_structure(self, base_structure, functions, cif_file):
        """
        Insert a structure with freeform functionalisation into the
        database. Will try to associate all the different functional
        groups contained.

        """

        new_mof_name = ".".join(functions)
        cif_file = "".join(cif_file)

        # Get an exisiting vesion of the structure if it exists
        structure = self.session.query(FreeFunctionalisedStructure).\
            filter(FreeFunctionalisedStructure.base_structure == base_structure).\
            filter(FreeFunctionalisedStructure.name == new_mof_name).\
            first()

        if structure is None:
            # Not in the database so initialise the structure and then add in
            # the functionalisations
            structure = FreeFunctionalisedStructure(base_structure,
                                                    new_mof_name,
                                                    cif_file)

            for unique_group in set(functions):
                group = self.groups.filter(FunctionalGroup.id == unique_group).first()
                structure.functional_groups.append(group)

            self.session.add(structure)
            # New structure inserted
        else:
            structure.cif_file = cif_file

        self.session.commit()

    def start_cif(self, structure_type, structure_id):
        """
        Return the cif file for the given structure and mark it as started
        in the database.

        """

        # bind the name as the extraction should be the same
        if structure_type == 'sym':
            debug("Looking for symmetry functionalised structure")
            Structure = SymFunctionalisedStructure
        elif structure_type == 'free':
            debug("Looking for free functionalised structure")
            Structure = FreeFunctionalisedStructure
        else:
            error("Unknown structure type %s" % structure_type)
            return None

        structure = self.session.query(Structure).\
            filter(Structure.id == structure_id).first()

        if structure is None:
            error("ID %i not found in database" % structure_id)
            return None
        else:
            # Mark it as started and send back the cif
            debug("Found structure %i" % structure_id)
            structure.status = STARTED
            self.session.commit()
            # TODO(tdaff): sqlite always returns a unicode object
            # test that this doesn't break across versions
            return structure.cif_file.encode('UTF-8')

    def store_results(self, structure_type, structure_id, structure):
        """
        Save the uptake data in an appropriate table. Will create tables if
        they do not exist.

        """

        # TODO(tdaff): more complete database to come
        if len(structure.guests) > 1:
            error("Database only deals with single guests at the moment")
            return

        # Just bind the guest we use (since there is only one here)
        guest = structure.guests[0]

        class GuestUptake(Base):
            """
            Guest uptake container, generated on the fly depending on guest
            and structure.

            """
            __tablename__ = '%s_%s_uptake' % (structure_type, guest.ident)

            # store all variations (necessary?)
            id = Column(Integer, primary_key=True)
            temperature = Column(Float)
            pressure = Column(Float)
            raw = Column(Float)
            raw_stdev = Column(Float)
            raw_supercell = Column(Integer)
            moluc = Column(Float)
            moluc_stdev = Column(Float)
            mmolg = Column(Float)
            mmolg_stdev = Column(Float)
            volvol = Column(Float)
            volvol_stdev = Column(Float)
            wtpc = Column(Float)
            wtpc_stdev = Column(Float)
            hoa = Column(Float)
            hoa_stdev = Column(Float)
            structure_id = Column(Integer, ForeignKey('%s_functionalised_structures.id' % structure_type))


        # create table for GuestUptake if doesn't already exist
        Base.metadata.create_all(self.engine)

        # Insert all the tp_points for the guest
        for tp_point in sorted(guest.uptake):
            # Instance row for each state point
            db_uptake = GuestUptake()
            db_uptake.structure_id = structure_id

            # Set the state point
            db_uptake.temperature = tp_point[0]
            db_uptake.pressure = tp_point[1][0]  # guest idx is 0 for one guest

            # keep the raw values just here
            # <N>, sd, supercell
            uptake = guest.uptake[tp_point]
            db_uptake.raw = uptake[0]
            db_uptake.raw_stdev = uptake[1]
            db_uptake.raw_supercell = uptake[2]

            # calculated values
            # normalise to unit cell
            uptake, stdev = (uptake[0]/uptake[2], uptake[1]/uptake[2])
            # molecules per unit cell
            db_uptake.moluc = uptake
            db_uptake.moluc_stdev = stdev
            # uptake in mmol/g
            db_uptake.mmolg = 1000*uptake/structure.weight
            db_uptake.mmolg_stdev = 1000*stdev/structure.weight
            # volumetric uptake
            db_uptake.volvol = (guest.molar_volume*uptake/(6.023E-4*structure.volume))
            db_uptake.volvol_stdev = (guest.molar_volume*stdev/(6.023E-4*structure.volume))
            # weight percent uptake
            db_uptake.wtpc = 100*(1 - structure.weight/(structure.weight + uptake*guest.weight))
            db_uptake.wtpc_stdev = 100*(1 - structure.weight/(structure.weight + stdev*guest.weight))
            # heat of adsorption
            hoa = guest.hoa[tp_point]
            db_uptake.hoa = hoa[0]
            db_uptake.hoa_stdev = hoa[1]

            # make sure to .commit() later
            self.session.add(db_uptake)

        # Tell the database the calculation is finished
        # bind the name as the extraction should be the same
        if structure_type == 'sym':
            Structure = SymFunctionalisedStructure
        elif structure_type == 'free':
            Structure = FreeFunctionalisedStructure
        else:
            error("Unknown structure type %s" % structure_type)
            return None

        structure = self.session.query(Structure).\
            filter(Structure.id == structure_id).first()

        if structure is None:
            error("ID %i not found in database" % structure_id)
            return None
        else:
            # Mark it as started and send back the cif
            debug("Set structure %i as finished in database" % structure_id)
            structure.status = FINISHED

        # finish up and save
        self.session.commit()
