#!/usr/bin/env python

"""
Sqlite backend to save generated structures in a database.

"""

from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, String, Text, ForeignKey
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()

# A secondary table is needed to link the structures to functionalisation sets
sym_associations = Table('sym_associations',
                         Base.metadata,
                         Column('functionalised_structure_id', Integer,
                                ForeignKey('sym_functionalised_structures.id')),
                         Column('functionalisation_id', Integer,
                                ForeignKey('functionalisations.id')))


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

    def __repr__(self):
        return "<SymFunctionalisedStructure('%s')>" % (self.fullname)

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

            print structure.fullname
            self.session.add(structure)
            # New structure inserted
        else:
            structure.cif_file = cif_file

        self.session.commit()
        print structure.id
