import datetime

import topology
import numpy as np
import pysheaf as ps
import pandas as pd
import datetime as dt
import json
import requests

class Sheaf():
    def __init__(self, open_set_list: topology.Topology):
        self.open_set_list = open_set_list

    def DistanceMeasurements(self, v1, v2):
        """Define how we measure distance between sheaf sections"""
        return np.linalg.norm(v2 - v1)

    def get_cell_name(self, open_set):
        cell_name = ""
        for x in open_set:
            cell_name = cell_name + x
        return cell_name

    def make_sheaf(self):
        print('Creating Sheaf')
        shf = ps.Sheaf()
        # Create sheaf cell for each open set in topology


        print('Creating Sheaf Cells')
        for i, open_set in enumerate(self.open_set_list):
            cell_name = self.get_cell_name(open_set)
            dimensionality = len(cell_name)
            dataType = f"{dimensionality}-Dimensional Vector"
            print(f'Creating sheaf cell {cell_name} of dimensionality {dimensionality}')
            shf.AddCell(cell_name, ps.Cell(dataType, self.DistanceMeasurements, dataDimension=dimensionality))

            #Set data for each open set
            shf.GetCell(cell_name).SetDataAssignment(
                ps.Assignment(dataType, topology.Topology().get_default_assignment(open_set)))


        #After all open sets are added to sheaf, add morphisms
        print('\n\n')
        print('***** Creating Restriction Morphisms *****')

        for i, open_set in enumerate(self.open_set_list):
            cell_name = self.get_cell_name(open_set)
            l = len(cell_name)
            domain_dimensionality = len(cell_name)
            domain_dataType = f"{domain_dimensionality}-Dimensional Vector"
            if l > 1:
                restriction_names_morphisms = []
                for x in range(0,l):
                    restriction = cell_name[0:x] + cell_name[x+1:l]
                    #restriction_names.append(restriction)
                    projection_morphism = np.eye(domain_dimensionality)
                    projection_morphism = np.delete(projection_morphism, x, 1)
                    #print(projection_morphism)
                    restriction_names_morphisms.append((restriction, projection_morphism))
                    #restriction_names = [cell_name[0:x] + cell_name[x+1:l] for x in range(0, l)]



                for (restriction, projection_morphism) in restriction_names_morphisms:
                    print(f' Creating restriction from {cell_name} to {restriction}')
                    codomain_dimensionality = len(restriction)
                    #projectionMorph = np.eye(domain_dimensionality, M=codomain_dimensionality)
                    print(projection_morphism)
                    codomain_dataType = f"{codomain_dimensionality}-Dimensional Vector"

                    shf.AddCoface(cell_name, restriction, \
                                  ps.Coface(domain_dataType, codomain_dataType, LinearMorphism(projection_morphism)))

        return shf

    def extend_sheaf(self, sheaf):
        for i, open_set in enumerate(self.open_set_list):
            cell_name = self.get_cell_name(open_set)
            sheaf.MaximallyExtendCell(cell_name)
        return sheaf

class SetMorphism():
  """A morphism in a subcategory of Set, described by a function object"""
  def __init__(self,fcn):
      self.fcn=fcn

  def __mul__(self,other): # Composition of morphisms
      return SetMorphism(lambda x : self.fcn(other.fcn(x)))

  def __call__(self,arg): # Calling the morphism on an element of the set
      return self.fcn(arg)

class LinearMorphism(SetMorphism):
    """A morphism in a category that has a matrix representation"""

    def __init__(self, matrix):
        SetMorphism.__init__(self, lambda x: np.dot(x, matrix))

    def __mul__(self, other):  # Composition of morphisms
        try:  # Try to multiply matrices.  This might fail if the other morphism isn't a LinearMorphism
            return LinearMorphism(np.dot(other.matrix, self.matrix))
        except AttributeError:
            return SetMorphism.__mul__(self, other)


class TimeSheaf:

    def __init__(self, df, max_depth=5, min_span=dt.timedelta(days=7)):

        self.df = df
        self.the_sheaf = ps.Sheaf()
        self.cell_names = []
        self.max_depth = max_depth  # maximum number of subdivisions of time period
        self.min_span = min_span  # smallest time span for grouping transactions

        self.time_topologize(min(df["Traded"]), max(df["Traded"]), 0)
        print(max(df["Traded"]) - min(df["Traded"]))

        for cn in self.cell_names:
            self.the_sheaf.MaximallyExtendCell(cn)



    def time_topologize(self, left_date, rite_date, depth):

        the_sheaf = self.the_sheaf

        #print('What is self.df.Transaction', self.df.Transaction == 'Purchase')
        #print(self.df.loc[(self.df.Transaction == 'Purchase')])
        the_transactions = self.df.loc[(self.df.Traded >= left_date)
                                       & (self.df.Traded <= rite_date) & (self.df.Transaction == 'Purchase')]
        print(the_transactions)
        if len(the_transactions) == 0:
            return None

        nm = str(left_date) + str(rite_date)
        self.cell_names.append(nm)

        the_sheaf.AddCell(nm, ps.Cell("Scalar"))
        the_sheaf.GetCell(nm).SetDataAssignment(ps.Assignment("Scalar",
                                                              self.get_data(the_transactions)))
        third = (rite_date - left_date) / 3
        if (depth < self.max_depth) and ((2 * third) > self.min_span):

            max_left = left_date + (third * 2)
            min_rite = left_date + third  # add from left (vs sub from rite) ensures overlap

            for L, R in [(left_date, max_left), (min_rite, rite_date)]:
                face = self.time_topologize(L, R, depth + 1)
                if face is not None:
                    #Morphism maps total growth over total range. Total range is leftmost thru rightmost date
                    the_sheaf.AddCoface(nm, face, ps.Coface("Scalar", "Scalar", lambda x: x))

        return nm

    def get_data(self, transactions):
        percent_changes = []
        for transaction in transactions:
            trade_date = transaction["Traded"]
            ticker = transaction['Ticker']

            initial_request_url = f"https://api.polygon.io/v1/open-close/{ticker}/{trade_date}?adjusted=true&apiKey=XrmofFN8e5WOGGBSx9iB0gkunwQN6R7d"
            json_data = json.loads(requests.get(initial_request_url).text)
            initial_price = json_data["close"]

            growth_period = datetime.timedelta(days=30)
            final_date = trade_date + growth_period
            final_request_url = f"https://api.polygon.io/v1/open-close/{ticker}/{final_date}?adjusted=true&apiKey=XrmofFN8e5WOGGBSx9iB0gkunwQN6R7d"
            json_data = json.loads(requests.get(final_request_url).text)
            final_price = json_data["close"]

            percent_change = (final_price-initial_price)/initial_price
            percent_changes.append(percent_change)
        # return transactions["excess_return"].dropna().mean()

        return percent_changes.mean()


if __name__ == "__main__":
    '''
    generators = topology.Topology().get_default_generator()
    osl = topology.Topology().generate_open_sets(generators)
    non_extended_sheaf = Sheaf(osl).make_sheaf()
    sheaf = Sheaf(osl).extend_sheaf(non_extended_sheaf)
    '''
    df = pd.read_excel("congress-trading-all.xlsx")
    df["Traded"] = pd.to_datetime(df["Traded"])
    '''
    for person in df["Name"].unique():
        print(person + " radius below")
        persondf = df.loc[(df.Name == person) & (df.Transaction == 'Purchase')]
        if len(persondf) > 0:
            if max(persondf["Traded"]) - min(persondf["Traded"]) >= dt.timedelta(days=7):
                t = TimeSheaf(df.loc[df.Name == person])

                print(t.the_sheaf.ComputeConsistencyRadius())
    '''
    """
    """
    t = TimeSheaf(df)

    for thres in [0., 0.1, 5]:
        print('Consistent stars at {} : {}'.format(thres,t.the_sheaf.ConsistentStarCollection(thres)))
    #print("Consistency Radius! " + str(sheaf.ComputeConsistencyRadius()))
