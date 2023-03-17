import numpy as np
from lmfit import Parameters, minimize
import pandas as pd

from models import gaussFit


def annulus_fit():


    df = pd.DataFrame({})