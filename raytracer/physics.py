"""physics.py"""
import numpy as np

def refract(direc=np.array([0,0,0]), normal=np.array([0,0,0]), n_1=0.0, n_2=0.0):
    """
    Implement Snell's law
    
    Args:
        direc: incident direction
        normal: surface normalised normal (away from n2 to n1)
        n_1: refractive index of 1
        n_2: refractive index of 2
        
    return:
        None: the ray is subject to TIR
        np.array(refracted_norm): normalised direction vector of refracted ray
    """
    direc_norm = direc / np.linalg.norm(direc)
    normal_norm = -(normal / np.linalg.norm(normal)) #minus sign for correct vector direction! (formula needs normal to point from n1 to n2)

    cos_theta_1 = np.dot(direc_norm, normal_norm) / \
        (np.linalg.norm(direc) * np.linalg.norm(normal))
    theta_1 = np.arccos(cos_theta_1)

    if np.sin(theta_1) > (n_2 / n_1):
        return None
    else:
        mu = n_1 / n_2
        refracted_norm = (mu * direc_norm) + \
                        normal_norm*np.sqrt(1-(mu**2)*(1-(np.dot(direc_norm, normal_norm))**2)) - \
                            (mu*normal_norm*np.dot(direc_norm, normal_norm))
        return np.array(refracted_norm)
