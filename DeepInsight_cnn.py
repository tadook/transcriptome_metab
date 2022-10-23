# Command Prompt
!git clone https://github.com/alok-ai-lab/pyDeepInsight.git
!pip install -r pyDeepInsight/requirements.txt

# Python
import sys
sys.path.append("pyDeepInsight/")
import pyDeepInsight.image_transformer

# create sample data

from sklearn.datasets import make_blobs
from matplotlib import pyplot
from pandas import DataFrame

# import host transcriptome dataset
import pandas as pd
trans_data = pd.read_csv('./normalized_counts.csv')
X = trans_data.drop(columns=['study_id','CPAPintubate','inpatient_hfo','severity','IntensiveTreatment'])
y = trans_data['severity']
print(X.shape, y.shape)

# generate a binary classification dataset (toy data)
# X, y = make_blobs(n_samples=400, centers=2, n_features=20000)
# print(X.shape, y.shape)

# scatter plot, dots colored by class value
# df = DataFrame(dict(x=X[:,3], y=X[:,4], label=y))
# colors = {0:'red', 1:'blue'}
# fig, ax = pyplot.subplots()
# grouped = df.groupby('label')

# for key, group in grouped:
#     group.plot(ax=ax, kind='scatter', x='x', y='y', label=key, color=colors[key])

# pyplot.show()

from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import RobustScaler

random_state=1515
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=23, stratify=y)

# mms = MinMaxScaler()
# X_train_norm = mms.fit_transform(X_train)
# X_test_norm = mms.transform(X_test)

# rs = RobustScaler()
# X_train_norm = rs.fit_transform(X_train)
# X_test_norm = rs.transform(X_test)

X_train_norm = X_train
X_test_norm = X_test

import numpy as np

le = LabelEncoder()
y_train_enc = le.fit_transform(y_train)
y_test_enc = le.transform(y_test)

le_mapping = dict(zip(le.transform(le.classes_), le.classes_))
num_classes = np.unique(y_train_enc).size

# import umap.umap_ as umap

# reducer = umap.UMAP(
#      n_components=2,
#      #min_dist=0.8,
#      metric='cosine',
#      n_jobs=-1
#  )

from sklearn.manifold import TSNE
distance_metric = 'cosine'

reducer = TSNE(
    n_components=2,
    metric=distance_metric,
    init='random',
    learning_rate='auto',
    n_jobs=-1
)

from pyDeepInsight import ImageTransformer

pixel_size = (224,224)
it = ImageTransformer(
    feature_extractor=reducer,
    pixels=pixel_size)
    
it.fit(X_train_norm, y=y_train, plot=True) # it takes 3 minites in case of "n:400, feature:20000"
X_train_img = it.transform(X_train_norm)
X_test_img = it.transform(X_test_norm)

X_train_img.shape, X_test_img.shape

import matplotlib.pyplot as plt

# fig, axs = plt.subplots(ncols=5, figsize=(3, 6))
# for i in range(5):
#   img = X_train_img[i]
#   label = y_train.iloc[i] # add ".iloc" which means integer based location.
#   print(f"{img.shape}, {img.min()}, {img.mean(axis=(0, 1))}, {img.max()}")
#   axs[i].imshow(img)
#   axs[i].set_title(f"sample {i}, label {label}, mean: {img.mean()}")

# plt.show()

# plot for each label

fix, axs = plt.subplots(ncols=5, nrows=num_classes, figsize=(40, 15))
for c in range(num_classes):
  indexes = np.where(y_train == c)[0]
  X_tmp = X_train_img[indexes]
  y_tmp = y_train.iloc[indexes]

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

output_img_dir = "outputs/input_imgs"

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

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
device

net = torch.hub.load('pytorch/vision:v0.10.0', 'squeezenet1_1', pretrained=True, verbose=False)
net.classifier[1] = nn.Conv2d(512, num_classes, kernel_size=(1,1), stride=(1,1))
net = net.to(device)

preprocess = transforms.Compose([
    transforms.ToTensor()
])

X_train_tensor = torch.stack([preprocess(img) for img in X_train_img]).float().to(device)
y_train_tensor = torch.from_numpy(le.fit_transform(y_train)).to(device)

X_test_tensor = torch.stack([preprocess(img) for img in X_test_img]).float().to(device)
y_test_tensor = torch.from_numpy(le.transform(y_test)).to(device)

batch_size = 64

trainset = TensorDataset(X_train_tensor, y_train_tensor)
trainloader = DataLoader(trainset, batch_size=batch_size, shuffle=True)

testset = TensorDataset(X_test_tensor, y_test_tensor)
testloader = DataLoader(testset, batch_size=batch_size, shuffle=False)

criterion = nn.CrossEntropyLoss()

# optimizer = optim.SGD(
#     net.parameters(),
#     lr=1e-04,
#     momentum=0.8,
#     weight_decay=1e-05
# )

optimizer = optim.Adam(net.parameters(), 1e-5)

def evaluate(net, dataloader):
  predicted_scores = np.empty(0)
  true_labels = np.empty(0)
  with torch.no_grad():
      net.eval()
      for i, data in enumerate(dataloader, 0):
          inputs, labels = data
          pred = torch.max(net(inputs),1)[1].cpu().detach().numpy()
          predicted_scores = np.append(predicted_scores, pred)
          true_labels = np.append(true_labels, labels.cpu().detach().numpy())

  metrics = {
      "accuracy": accuracy_score(true_labels, predicted_scores),
      "auroc": roc_auc_score(true_labels, predicted_scores)
  }

  return predicted_scores, true_labels, metrics

# training loop
from sklearn.metrics import roc_auc_score

num_epochs = 700
for epoch in range(num_epochs):
    net.train()

    running_loss = 0.0
    for i, data in enumerate(trainloader, 0):
        # get the inputs; data is a list of [inputs, labels]
        inputs, labels = data

        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        outputs = net(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        running_loss += loss.item()

    # print epoch statistics
    if not (epoch % 10):
        print(f'[Epoch {epoch}] loss: {running_loss / len(X_train_tensor) * batch_size:.3f}')

        _, _, train_metrics = evaluate(net, trainloader)
        print("train", train_metrics)

        _, _, test_metrics = evaluate(net, testloader)
        print("test", test_metrics)
  
print(f'[{epoch}] loss: {running_loss / len(X_train_tensor) * batch_size:.3f}')

# run predictions on training samples
train_predicted, train_true, train_metrics = evaluate(net, trainloader)
print(train_metrics)

# run predictions on test samples
test_predicted, test_true, test_metrics = evaluate(net, testloader)
print(test_metrics)

# # run predictions on training samples
# train_predicted = np.empty(0)
# train_true = np.empty(0)
# with torch.no_grad():
#     net.eval()
#     for i, data in enumerate(trainloader, 0):
#         inputs, labels = data
#         pred = torch.max(net(inputs),1)[1].cpu().detach().numpy()
#         train_predicted = np.append(train_predicted, pred)
#         train_true = np.append(train_true, labels.cpu().detach().numpy())

# # run predictions on test samples
# test_predicted = np.empty(0)
# test_true = np.empty(0)
# with torch.no_grad():
#     net.eval()
#     for i, data in enumerate(testloader, 0):
#         inputs, labels = data
#         pred = torch.max(net(inputs),1)[1].cpu().detach().numpy()
#         test_predicted = np.append(test_predicted, pred)
#         test_true = np.append(test_true, labels.cpu().detach().numpy())

print(f"The train accuracy was {roc_auc_score(train_true, train_predicted):.3f}")
print(f"The test accuracy was {roc_auc_score(test_true, test_predicted):.3f}")



# Grad-CAM

!pip install grad-cam

from pytorch_grad_cam import GradCAM, HiResCAM, ScoreCAM, GradCAMPlusPlus, AblationCAM, XGradCAM, EigenCAM, FullGrad
from pytorch_grad_cam.utils.model_targets import ClassifierOutputTarget
from pytorch_grad_cam.utils.image import show_cam_on_image
from torchvision.models import resnet50

model = resnet50(pretrained=True)
target_layers = [model.layer4[-1]]
input_tensor = # Create an input tensor image for your model..
# Note: input_tensor can be a batch tensor with several images!

# Construct the CAM object once, and then re-use it on many images:
cam = GradCAM(model=model, target_layers=target_layers, use_cuda=args.use_cuda)

# You can also use it within a with statement, to make sure it is freed,
# In case you need to re-create it inside an outer loop:
# with GradCAM(model=model, target_layers=target_layers, use_cuda=args.use_cuda) as cam:
#   ...

# We have to specify the target we want to generate
# the Class Activation Maps for.
# If targets is None, the highest scoring category
# will be used for every image in the batch.
# Here we use ClassifierOutputTarget, but you can define your own custom targets
# That are, for example, combinations of categories, or specific outputs in a non standard model.

targets = [ClassifierOutputTarget(281)]

# You can also pass aug_smooth=True and eigen_smooth=True, to apply smoothing.
grayscale_cam = cam(input_tensor=input_tensor, targets=targets)

# In this example grayscale_cam has only one image in the batch:
grayscale_cam = grayscale_cam[0, :]
visualization = show_cam_on_image(rgb_img, grayscale_cam, use_rgb=True)