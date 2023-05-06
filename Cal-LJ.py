import MDAnalysis as mda
import numpy as np
import os
import matplotlib.pyplot as plt

# 定义MXene膜的原子名称和LJ参数
mxene_atoms = {
    'Ti1': {'sigma': 2.829, 'epsilon': 0.07108},
    'Ti2': {'sigma': 2.829, 'epsilon': 0.07108},
    'C1': {'sigma': 3.401, 'epsilon': 0.43897},
    'C2': {'sigma': 3.401, 'epsilon': 0.43897}
}

# 定义要计算相互作用能的原子和对应的LJ参数
atom_lj_params = {
    'Cn': {'sigma': 3.730, 'epsilon': 1.22964},
    #'O': {'sigma': 0.31, 'epsilon': 1.65}
}

# 定义LB混合规则
def lb_mixing_rule(sigma1, sigma2, epsilon1, epsilon2):
    sigma = (sigma1 + sigma2) / 2
    epsilon = np.sqrt(epsilon1 * epsilon2)
    return sigma, epsilon

# 创建输出文件
output_file = open("LJ.txt", "w")


for i in range(135, 171):
    # 生成文件名
    gro_file = f"model/CH4-{i/10:.1f}.gro"
    
    # 读取gro文件
    u = mda.Universe(gro_file)
    
    # 计算每种原子对的相互作用能并求和
    total_energy = 0.0
    for atom1 in atom_lj_params:
        for atom2 in mxene_atoms:
            # 获取原子的LJ参数
            lj_params1 = atom_lj_params[atom1]
            lj_params2 = mxene_atoms[atom2]
            
            # 计算不同原子对的sigma和epsilon参数
            sigma, epsilon = lb_mixing_rule(lj_params1['sigma'], lj_params2['sigma'],
                                             lj_params1['epsilon'], lj_params2['epsilon'])
            
            # 获取原子的索引
            atom1_indices = u.select_atoms('name ' + atom1).indices
            atom2_indices = u.select_atoms('name ' + atom2).indices
            
            # 计算相互作用能
            positions = u.atoms.positions
            r = positions[atom2_indices, np.newaxis, :] - positions[atom1_indices]
            distances = np.sqrt(np.sum(r**2, axis=-1))
            total_energy += np.sum(4 * epsilon * ((sigma / distances)**12 - (sigma / distances)**6))
    
    # 输出结果到文件
    #output_file.write(f"Interaction energies for {gro_file}:\n")
    #output_file.write(f"Total LJ energy = {total_energy}\n\n")
    output_file.write("{:.2f} {:.6f}\n".format(i/10, total_energy))
    #output_file.write(f"{total_energy}\n")

# 关闭输出文件
output_file.close()

# 画图
data = np.loadtxt("LJ.txt")
x = data[:, 0]
y = data[:, 1]
plt.plot(x, y, 'o-', color='red')
plt.show()
