from knoten import csm, utils
import spiceypy as spice
import numpy as np

class Illuminator:

    def __init__(self, position=None, velocity=None):
        self.position = position
        self.velocity = velocity

    def get_position_from_csm_sensor(self, sensor, ground_pt):
        sunEcefVec = sensor.getIlluminationDirection(ground_pt)
        self.position = utils.Point(ground_pt.x - sunEcefVec.x, ground_pt.y - sunEcefVec.y, ground_pt.z - sunEcefVec.z)
        return self.position 
    
    # def get_position_from_time(self, sensor_state):
    #     return 0
    
    # def get_velocity(self, sensor_state):
    #     return 0


        