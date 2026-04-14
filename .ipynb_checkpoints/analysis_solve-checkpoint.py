import pandas as pd
import numpy as np

class AnalysisFunc:
    def __init__(self):
        pass
    
    def situation_1(self, x, M, E, I, l):
        """情况1: 纯弯矩作用下的挠度"""
        return -M * x * x / 2 / E / I
    
    def situation_2(self, x, F, E, I, l):
        """情况2: 自由端受集中力作用的悬臂梁"""
        return -F * x * x / 6 / E / I * (3 * l - x)
    
    def situation_3(self, x, F, a, E, I, l):
        """情况3: 悬臂梁在距固定端a处受集中力作用"""
        if x <= a:
            return -F * x * x / 6 / E / I * (3 * a - x)
        else:
            return -F * a * a / 6 / E / I * (3 * x - a)
    
    def situation_4(self, x, q, E, I, l):
        """情况4: 均布载荷作用下的挠度"""
        return -q * x * x / 24 / E / I * (x * x - 4 * l * x + 6 * l * l)
    
    def generate_situation_1_data(self, M, E, I, l, q, F, a, num_points=2000):
        """生成situation_1的数据点"""
        x_points = np.linspace(0, l, num_points)
        y_points = [self.situation_1(x, M, E, I, l) for x in x_points]
        
        df = pd.DataFrame({
            'x': x_points,
            'w': y_points
        })
        return df
    
    def generate_situation_2_data(self, M, E, I, l, q, F, a, num_points=2000):
        """生成situation_2的数据点"""
        x_points = np.linspace(0, l, num_points)
        y_points = [self.situation_2(x, F, E, I, l) for x in x_points]
        
        df = pd.DataFrame({
            'x': x_points,
            'w': y_points
        })
        return df
    
    def generate_situation_3_data(self, M, E, I, l, q, F, a, num_points=2000):
        """生成situation_3的数据点"""
        x_points = np.linspace(0, l, num_points)
        y_points = [self.situation_3(x, F, a, E, I, l) for x in x_points]
        
        df = pd.DataFrame({
            'x': x_points,
            'w': y_points
        })
        return df
    
    def generate_situation_4_data(self, M, E, I, l, q, F, a, num_points=2000):
        """生成situation_4的数据点"""
        x_points = np.linspace(0, l, num_points)
        y_points = [self.situation_4(x, q, E, I, l) for x in x_points]
        
        df = pd.DataFrame({
            'x': x_points,
            'w': y_points
        })
        return df
