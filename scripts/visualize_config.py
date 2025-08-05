import os
import re
import sys
import matplotlib.pyplot as plt
from tabulate import tabulate
from mpl_toolkits.mplot3d import Axes3D

def parse_log_file(filepath):
    config_data = {}
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Extract key-value pairs for configuration
            if "Bank Organization" in line:
                config_data["Bank Organization"] = re.findall(r'\d+ x \d+', line)[0]
            elif "Mat Organization" in line:
                config_data["Mat Organization"] = re.findall(r'\d+ x \d+', line)[0]
            elif "Subarray Size" in line:
                config_data["Subarray Size"] = re.findall(r'\d+ Rows x \d+ Columns', line)[0]
            elif "Senseamp Mux" in line:
                config_data["Senseamp Mux"] = re.findall(r'\d+', line)[0]
            elif "Output Level-1 Mux" in line:
                config_data["Output Level-1 Mux"] = re.findall(r'\d+', line)[0]
            elif "Output Level-2 Mux" in line:
                config_data["Output Level-2 Mux"] = re.findall(r'\d+', line)[0]
            elif "Total Area" in line:
                config_data["Total Area"] = re.findall(r'[\d.]+um x [\d.]+um', line)[0]
            elif "Area Efficiency" in line:
                config_data["Area Efficiency"] = re.findall(r'[\d.]+%', line)[0]
            elif "Read Latency" in line:
                config_data["Read Latency"] = re.findall(r'[\d.]+ns', line)[0]
            elif "Write Latency" in line:
                config_data["Write Latency"] = re.findall(r'[\d.]+ns', line)[0]
            elif "Read Bandwidth" in line:
                config_data["Read Bandwidth"] = re.findall(r'[\d.]+GB/s', line)[0]
            elif "Write Bandwidth" in line:
                config_data["Write Bandwidth"] = re.findall(r'[\d.]+MB/s', line)[0]
    return config_data

def draw_layout(config_data):
    # Parse organization sizes
    bank_rows, bank_cols = map(int, config_data["Bank Organization"].split(' x '))
    mat_rows, mat_cols = map(int, config_data["Mat Organization"].split(' x '))
    sub_rows, sub_cols = map(int, re.findall(r'(\d+) Rows x (\d+) Columns', config_data["Subarray Size"])[0])

    fig, ax = plt.subplots(figsize=(10, 8))
    mat_size = 1.0
    bank_spacing = mat_cols * mat_size + 1.0
    mat_spacing = mat_size + 0.2

    for b_row in range(bank_rows):
        for b_col in range(bank_cols):
            bank_x = b_col * bank_spacing
            bank_y = b_row * bank_spacing
            # Calculate bank center for routing
            bank_center_x = bank_x + (mat_cols * mat_spacing) / 2
            bank_center_y = bank_y + (mat_rows * mat_spacing) / 2
            # Draw mats inside each bank
            for m_row in range(mat_rows):
                for m_col in range(mat_cols):
                    mat_x = bank_x + m_col * mat_spacing
                    mat_y = bank_y + m_row * mat_spacing
                    rect = plt.Rectangle((mat_x, mat_y), mat_size, mat_size, fill=True, color='lightblue', edgecolor='blue')
                    ax.add_patch(rect)
                    ax.text(mat_x+mat_size/2, mat_y+mat_size/2, f'M', ha='center', va='center', fontsize=8)
                    # Draw routing from mat center to bank center
                    ax.plot([mat_x+mat_size/2, bank_center_x], [mat_y+mat_size/2, bank_center_y], color='red', lw=1, linestyle='--')
            # Draw bank rectangle
            bank_rect = plt.Rectangle((bank_x-0.2, bank_y-0.2), mat_cols*mat_spacing+0.4, mat_rows*mat_spacing+0.4, fill=False, edgecolor='black', lw=2)
            ax.add_patch(bank_rect)
            ax.text(bank_x + (mat_cols*mat_spacing)/2, bank_y + (mat_rows*mat_spacing)+0.3, f'Bank', ha='center', va='bottom', fontsize=10, color='black')

    ax.set_xlim(0, bank_cols * bank_spacing)
    ax.set_ylim(0, bank_rows * bank_spacing + 1)
    ax.set_aspect('equal')
    ax.set_title("Memory Architecture Layout\n(Banks > Mats > Routing)")
    ax.axis('off')
    plt.show()

def get_num_layers(filepath):
    # Try to infer layers from config filename or log content
    with open(filepath, 'r') as file:
        for line in file:
            if "sample_3DReRAM.cfg" in line:
                return 4  # Example: set to 4 layers for 3D config
            # You can add more robust parsing here if layers info is available
    return 1  # Default to 1 layer if not found

def draw_layout_3d(config_data, num_layers=1):
    bank_rows, bank_cols = map(int, config_data["Bank Organization"].split(' x '))
    mat_rows, mat_cols = map(int, config_data["Mat Organization"].split(' x '))
    sub_rows, sub_cols = map(int, re.findall(r'(\d+) Rows x (\d+) Columns', config_data["Subarray Size"])[0])

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    mat_size = 1.0
    bank_spacing = mat_cols * mat_size + 1.0
    mat_spacing = mat_size + 0.2
    layer_spacing = 2.0

    for layer in range(num_layers):
        z = layer * layer_spacing
        for b_row in range(bank_rows):
            for b_col in range(bank_cols):
                bank_x = b_col * bank_spacing
                bank_y = b_row * bank_spacing
                bank_center_x = bank_x + (mat_cols * mat_spacing) / 2
                bank_center_y = bank_y + (mat_rows * mat_spacing) / 2
                bank_center_z = z
                for m_row in range(mat_rows):
                    for m_col in range(mat_cols):
                        mat_x = bank_x + m_col * mat_spacing
                        mat_y = bank_y + m_row * mat_spacing
                        ax.bar3d(mat_x, mat_y, z, mat_size, mat_size, 0.5, color='lightblue', edgecolor='blue', alpha=0.7)
                        ax.text(mat_x+mat_size/2, mat_y+mat_size/2, z+0.3, f'M', ha='center', va='center', fontsize=7)
                        # Routing
                        ax.plot([mat_x+mat_size/2, bank_center_x], [mat_y+mat_size/2, bank_center_y], [z+0.25, bank_center_z+0.25], color='red', lw=0.5, linestyle='--')
                # Bank label
                ax.text(bank_center_x, bank_center_y, z+0.7, f'Bank L{layer+1}', ha='center', va='bottom', fontsize=9, color='black')

    ax.set_title("3D Memory Architecture Layout\n(Banks > Mats > Routing > Layers)")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Layer')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python visualize_config.py <log_path>")
        sys.exit(1)
    log_path = sys.argv[1]
    if not os.path.exists(log_path):
        print(f"Log file not found at {log_path}")
    else:
        config = parse_log_file(log_path)
        num_layers = get_num_layers(log_path)
        if num_layers > 1:
            draw_layout_3d(config, num_layers)
        else:
            draw_layout(config)
