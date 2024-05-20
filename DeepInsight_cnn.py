# Command Prompt
# cd /Users/tao/Library/CloudStorage/Dropbox/Tao/1_HasegawaLab/papers/transcriptome_metabolome/work/Github
!git clone https://github.com/alok-ai-lab/pyDeepInsight.git
!pip install -r ../pyDeepInsight/requirements.txt
!pip install seaborn
!pip install imblearn
!pip install plotly_express==0.4.0
!pip install rpy2

# Python
import sys
sys.path.append("../pyDeepInsight/")
import pyDeepInsight.image_transformer

from sklearn.datasets import make_blobs
from matplotlib import pyplot
from pandas import DataFrame

# import dataset
import pandas as pd
tra_met_data = pd.read_csv('../trans_metab_int_rev0506.csv') 
X = tra_met_data.drop(columns=['study_id','CPAPintubate','inpatient_hfo','severity','IntensiveTreatment','intake_sex','Age_mo'])
y = tra_met_data['severity']
genes = tra_met_data.iloc[:, 7:5558].columns.to_numpy() # 7:6706 for sensitivity analysis
print(X.shape, y.shape) 

# make train/test data
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import RobustScaler
from pyDeepInsight import Norm2Scaler

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=24, stratify=y)

# save X_train/test and y_train/test for check with other models
# X_train.to_csv("outputs/X_train_412.csv")
# y_train.to_csv("outputs/y_train_412.csv")
# X_test.to_csv("outputs/X_test_412.csv")
# y_test.to_csv("outputs/y_test_412.csv")

# SMOTE
from imblearn.over_sampling import SMOTE

sm = SMOTE(random_state=42)
X_train, y_train = sm.fit_resample(X_train, y_train)

X_train_norm = X_train
X_test_norm = X_test

X_train_norm = X_train_norm.fillna(0)
X_test_norm = X_test_norm.fillna(0)

# Label encode and Mapping
import numpy as np

le = LabelEncoder()
y_train_enc = le.fit_transform(y_train)
y_test_enc = le.transform(y_test)

le_mapping = dict(zip(le.transform(le.classes_), le.classes_))
num_classes = np.unique(y_train_enc).size

# tSNE

from sklearn.manifold import TSNE
distance_metric = 'cosine'

reducer = TSNE(
    n_components=2,
    metric=distance_metric,
    init='random',
    learning_rate='auto',
    n_jobs=-1,
    verbose=True,
    random_state=42
)

# Image Transform
from pyDeepInsight import ImageTransformer

pixel_size = (224, 224) 
it = ImageTransformer(
    feature_extractor=reducer,
    pixels=pixel_size)

it.fit(X_train, y=y_train, plot=True) 

X_train_img = it.transform(X_train_norm)
X_test_img = it.transform(X_test_norm)

X_train_img.shape, X_test_img.shape

pd.to_pickle(it,'models/image/it_0506.pkl')


# Convert images to ones with colors

coords = it.coords()

## train data
X_img = X_train_img
new_img = np.zeros_like(X_img)

for k in range(0, 5550): # All variables (6698 for sensitivity analysis)
    x = coords[k][0]
    y = coords[k][1]

    if k <= 127: # Red: metabolite 128
        new_img[:, x, y, 0] = X_img[:, x, y, 0]

    elif k >= 128 and k <= 715: # Green: transcriptome 588
        new_img[:, x, y, 1] = X_img[:, x, y, 1]

    else: # Blue: metabolite * transcriptome 
        new_img[:, x, y, 2] = X_img[:, x, y, 2]

X_train_img = new_img

## test data
X_img = X_test_img
new_img = np.zeros_like(X_img)

for k in range(0, 5550): # All variables (6698 for sensitivity analysis)
    x = coords[k][0]
    y = coords[k][1]

    if k <= 127: # Red: metabolite 128
        new_img[:, x, y, 0] = X_img[:, x, y, 0]

    elif k >= 128 and k <= 715: # Green: transcriptome 588
        new_img[:, x, y, 1] = X_img[:, x, y, 1]

    else: # Blue: metabolite * transcriptome
        new_img[:, x, y, 2] = X_img[:, x, y, 2]

X_test_img = new_img

## control image data
new_img = np.zeros_like(X_img)
new_img2 = np.ones_like(X_img)

for k in range(0, 5550): # All variables (6698 for sensitivity analysis)
    x = coords[k][0]
    y = coords[k][1]

    if k <= 127: # Red: metabolite 128
        new_img[:, x, y, 0] = new_img2[:, x, y, 0]

    elif k >= 128 and k <= 715: # Green: transcriptome 588
        new_img[:, x, y, 1] = new_img2[:, x, y, 1]

    else: # Blue: metabolite * transcriptome
        new_img[:, x, y, 2] = new_img2[:, x, y, 2]

X_cont_img = new_img

# save X_train/test_img and y_train/test for Google Colab
import numpy
numpy.save('outputs/tramet_tsne_0506/X_train_img.npy', X_train_img)
numpy.save('outputs/tramet_tsne_0506/X_test_img.npy', X_test_img)
y_train.to_pickle('outputs/tramet_tsne_0506/y_train.pkl')
y_test.to_pickle('outputs/tramet_tsne_0506/y_test.pkl')

# View overall feature overlap
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

fdm = it.feature_density_matrix()
fdm[fdm == 0] = np.nan

plt.figure(figsize=(10, 7.5))

ax = sns.heatmap(fdm, cmap="viridis", linewidths=0., 
                 linecolor="lightgrey", square=True)
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
for _, spine in ax.spines.items():
    spine.set_visible(True)
_ = plt.title("Genes per pixel")

plt.show()

# plot for each label
fix, axs = plt.subplots(ncols=5, nrows=num_classes, figsize=(40, 15))
for c in range(num_classes):
  indexes = np.where(y_train == c)[0]
  X_tmp = X_train_img[indexes]
  y_tmp = y_train.iloc[indexes]

  for i in range(5):
    img = X_tmp[i+15]
    label = y_tmp.iloc[i+15]
    print(f"{img.shape}, {img.min()}, {img.mean(axis=(0, 1))}, {img.max()}")
    axs[c, i].imshow(img)
    axs[c, i].set_title(f"sample {i}, label {label}")

fix.tight_layout() 

plt.show()

# mean image for each label (you have to use y_train data before SMOTE)
import numpy as np
import matplotlib.pyplot as plt

indexes_train_1 = np.where(y_train == 1)[0]
indexes_train_0 = np.where(y_train == 0)[0]
indexes_test_1 = np.where(y_test == 1)[0]
indexes_test_0 = np.where(y_test == 0)[0]

X_tmp_1 = np.concatenate([X_train_img[indexes_train_1], X_test_img[indexes_test_1]], axis=0)
X_tmp_0 = np.concatenate([X_train_img[indexes_train_0], X_test_img[indexes_test_0]], axis=0)

mean_img_1 = np.mean(X_tmp_1, axis=0)
mean_img_0 = np.mean(X_tmp_0, axis=0)

fig, axs = plt.subplots(ncols=2, figsize=(10, 5))

axs[0].imshow(mean_img_1)
axs[0].set_title("Mean image for label 1")

axs[1].imshow(mean_img_0)
axs[1].set_title("Mean image for label 0")

plt.show()


# plot the control image
fix, axs = plt.subplots(ncols=5, nrows=num_classes, figsize=(40, 15))
for c in range(num_classes):
  indexes = np.where(y_test == c)[0]
  X_tmp = X_cont_img[indexes]
  y_tmp = y_test.iloc[indexes]

  for i in range(5):
    img = X_tmp[i]
    label = y_tmp.iloc[i]
    print(f"{img.shape}, {img.min()}, {img.mean(axis=(0, 1))}, {img.max()}")
    axs[c, i].imshow(img)
    axs[c, i].set_title(f"sample {i}, label {label}")

fix.tight_layout() 

plt.show()

# save images
import os
from PIL import Image

output_img_dir = "outputs/imgs_tramet_0506_82"

for split, X, y in zip(['train', 'test'], [X_train_img, X_test_img], [y_train, y_test]):
  for c in range(num_classes):
    indexes = np.where(y == c)[0]
    X_tmp = X[indexes]
    y_tmp = y.iloc[indexes]  

    output_img_class_dir = os.path.join(output_img_dir, split, str(c))
    os.makedirs(output_img_class_dir, exist_ok=True)  

    print(f"saving {len(y_tmp)} {split} images for class {c} to {output_img_class_dir}")  

    for i in range(len(y_tmp)):
      idx = indexes[i]
      img = X_tmp[i]
      label = y_tmp.iloc[i]

      img_path = os.path.join(output_img_class_dir, f"idx{idx}_class{label}.png")
      Image.fromarray((img * 255).astype(np.uint8)).save(img_path)


# CNN 
import torch
import torchvision.transforms as transforms
from torch.utils.data import TensorDataset, DataLoader
import torch.nn as nn
import torch.optim as optim

from sklearn.metrics import accuracy_score

import warnings; 
warnings.simplefilter('ignore')

# load the trained model
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
net = torch.hub.load('pytorch/vision:v0.10.0', 'squeezenet1_1', pretrained=True, verbose=False)
num_classes = 2 # binominal 0/1 -> 2
net.classifier[1] = nn.Conv2d(512, num_classes, kernel_size=(1,1), stride=(1,1))
net = net.to(device)
net.load_state_dict(torch.load("models/cnn/last_model.pth", map_location=torch.device('cpu'))) # in GPU environment, we don't need map_location option

# load X_train/test_img and y_train/test
import pickle
import numpy as np
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
X_train_img = np.load('outputs/tramet_tsne_0412/X_train_img.npy')
X_test_img = np.load('outputs/tramet_tsne_0412/X_test_img.npy')
with open('outputs/tramet_tsne_0412/y_train.pkl', "rb") as fh:
  y_train = pickle.load(fh)
with open('outputs/tramet_tsne_0412/y_test.pkl', "rb") as fh:
  y_test = pickle.load(fh)

preprocess = transforms.Compose([
    transforms.ToTensor()
])

le = LabelEncoder()
y_train_enc = le.fit_transform(y_train)
y_test_enc = le.transform(y_test)
le_mapping = dict(zip(le.transform(le.classes_), le.classes_))

X_train_tensor = torch.stack([preprocess(img) for img in X_train_img]).float().to(device)
y_train_tensor = torch.from_numpy(le.fit_transform(y_train)).to(device)

X_test_tensor = torch.stack([preprocess(img) for img in X_test_img]).float().to(device)
y_test_tensor = torch.from_numpy(le.transform(y_test)).to(device)

# load the image
import pickle
with open('models/image/it_0412.pkl', "rb") as fh:
  it = pickle.load(fh)

# Deep Feature: CAM-based feature selection
from pyDeepInsight import CAMFeatureSelector

## Step 1 - CAMFeatureSelector object
cm_method='GradCAM'
camfs = CAMFeatureSelector(
    model=net,
    it=it,
    cam_method=cm_method
)

## Step 2 - Compute Class-Specific CAMs
fl_method = "max"
class_cam = camfs.calculate_class_activations(X_train_tensor, y_train_tensor, batch_size=100, flatten_method=fl_method)

## Step 3 - Select Class-Specific Features

## Added modification to .select_class_features (Line191: pyDeepInsight/pyDeepInsight/feature_selection.py)
# def select_class_features(self, cams: Dict[int, np.ndarray],
#                             threshold: float = 0.6) -> Dict[int, np.ndarray]:
#     class_idx = {}
#     for cat, cam in cams.items():
#         cam_pass = np.stack(np.where(cam >= threshold)).T
#         it_pass = np.where(
#             (self.feature_coords == cam_pass[:, None]).all(-1) #### Modified ####
#         )[1]
#         class_idx[cat] = it_pass
#     return class_idx
                 
fs_threshold = 0.001
feat_idx = camfs.select_class_features(cams=class_cam, threshold=fs_threshold)
feat_idx

from pytorch_grad_cam.utils.image import show_cam_on_image
from matplotlib import pyplot as plt

def cam_image(X, y, cam, fs, threshold):
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(2, 1),
                            constrained_layout=True, squeeze=False) # added "squeeze=False" for binomial outcome
    for cat in np.unique(y):
        row = cat // 2 # cat // 4
        col = cat # cat % 4
        cat_idx = np.where(y == cat)[0]
        X_cat = X[cat_idx,:,:,:].detach().mean(dim=0).cpu().numpy()
        cam_cat = cam[cat].copy()
        cam_cat[cam_cat <= threshold] = 0
        visualization = show_cam_on_image(
            np.transpose(X_cat, (1,2,0)),
            cam_cat,
            use_rgb=True
        )
        _ = axs[row, col].imshow(visualization)
        axs[row, col].text(0,0,le_mapping[cat],c="white",ha="left",va="top",weight="bold",size="x-large")
        axs[row, col].text(224,224,f"{fs[cat].shape[0]} genes",c="white",ha="right",va="bottom",weight="bold",size="large")
        axs[row, col].axis('off')
    return fig, axs

_ = cam_image(X_train_tensor, y_train_tensor.detach().cpu().numpy(), class_cam, feat_idx, fs_threshold)

plt.show()

## Step 4 - Extract Feature Names

### Category
for cat, idx in feat_idx.items():
    feature_names = genes[idx] 
    print(f"{idx.shape[0]:5} features selected for {le_mapping[cat]:4}: {', '.join(feature_names[1:10])}...")
    
    feat = pd.DataFrame()

feature_severe_0 = genes[list(feat_idx.items())[0][1]] 
feature_severe_1 = genes[list(feat_idx.items())[1][1]] 

feature_severe_0
feature_severe_1

# Step 5 - Extract the importance of each feature by fetching maximum CAM values across all images

# Initialize a dictionary to store the importance of each feature
feature_importance = {gene: 0 for gene in genes}

# Define the range of thresholds to consider
thresholds = np.arange(0.001, 1.001, 0.001)

 # Select the features that exceed the threshold
for threshold in thresholds:  
    feat_idx = camfs.select_class_features(cams=class_cam, threshold=threshold)
    for cat, idx in feat_idx.items():
        for feature in idx:
            # Update the importance of the feature
            feature_importance[genes[feature]] = threshold

# Convert the dictionary to a DataFrame
feature_importance_df = pd.DataFrame.from_dict(feature_importance, orient='index', columns=['importance'])

# Save the DataFrame to a CSV file
feature_importance_df.to_csv('outputs/feature_importance.csv')


### Table
feat = pd.DataFrame()
for cat, idx in feat_idx.items():
    feature_names = genes[idx]
    feat = pd.concat([feat, pd.DataFrame({'severity':le_mapping[cat], 'gene':feature_names})])
fdf = feat.assign(selected=1).pivot(index='severity', columns='gene', values="selected").fillna(0).astype(int)

pd.DataFrame(
    np.matmul(fdf.values,fdf.T.values),
    index=fdf.index.values,
    columns=fdf.index.values
)

#### save the feature
import csv

f = open('../feature_severe.csv', 'w',newline="")
writer = csv.writer(f)
for w in range(len(feature_severe_0)):
 writer.writerow([feature_severe_0[w]])
f.close()
