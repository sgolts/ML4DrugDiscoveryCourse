# Train deep neural networks

In this lab you will develop and train a deep-neural network to
predict docking scores from molecular fingerprints


# Learning objectives

* Familiarity with basic Pytorch deep-learning models and training
* How to approach fitting hyperparameters and analyzing performance

# Steps
This lab will require using a GPU to run train the deep learning models, so I recommend
running this on the compute cluster. You could do the initial model development locally
or in a non-GPU node and then debug and run the model interactively or in batch mode
using a GPU node.


### 1 Setup Weights and Biases to do experiment tracking and set up your environment
[Weights and biases](https://wandb.ai) is a website/service for experiment tracking.

Follow the quickstart instructions

  1) Create an account and login to [https://wandb.ai](https://wandb.ai), create an API key, and create a project e.g.
     `BIOINF595w26_lab_08`.


  2) In a conda environment in your compute environment, install the python packages

    pip install wandb

  3) Login on the command line and when prompted give the API key.

    wandb login

  4) Load modules and install deep-learning packages. Consider loading the cuda and cudnn modules

    module load cuda cudnn

  Install pytorch

    pip3 install torch

  On a compute node with a GPU check that it you are able to access the GPU

    nvidia-smi
    python -c "import torch; print(f'Is torch available? {torch.cuda.is_available()}')"

### 2 Download and split one of the LSD Docking datasets

   1) Similar to the data we used for the supervised ML data in lab 04, gather data from one of the
      `IrwinLab` HuggingFace datasets, and sample a small and medium datasets with `10_000` and `100_000`
      compounds each.

   2) Split the data e.g. into e.g. 60% for train, 20% for validation, and 20% for test using
      `sklearn.module_selection.train_test_split`.


### 3 Create FingerprintNN Network

Create a python file `model.py` with the following parts

  1) Implement the `FingerprintNN` class that derives from `torch.nn.Module`. This class defines
     the pytorch model and stores the model parameters. It supports initialization through
     the `reset_parameters` and evaluation through the `forward` function. It should implement
     the following functions:

     a) Define the `__init__` function to take parameters

         * input_dim: integer dimension of the input features
         * hidden_dim: integer dimension of the hidden layers
         * n_layers: integer number of layers in the network

        and define a `self.model` to be a `nn.Sequential(...)` object which takes in a list
        of modules that are evaluated sequentially. Initialize it with
        `nn.Linear(...)` followed by `nn.ReLU()` repeating `n_layers` times. The dimensions should be
        `input_dim`, `hidden_dim`, ..., `hidden_dim`, `1`.

     b) Define the `forward` function to take in `x`, which should be a `torch.tensor` with dimensions
        `[batch_size, input_dim]` and applies `self.model` as a function to `x`. Next remove dimensions
        that are dimension 1 by calling the `.sqeeze()` function on the resulting object
        to return it as a `torch.tensor` of dimension `[batch_size]`.

     c) Define the `reset_parameters` function that initializes the parameters like this:

        for m in self.modules():
            if isinstance(m, nn.Linear):
                 nn.init.xavier_uniform_(m.weight)
                 nn.init.zeros_(m.bias)

  2) Implement the `FingerprintDataset` class that derives from
     `torch.utils.data.Dataset` and can supply the model with data for
     training, validation, and testing. Here we're going to take in a
     Pandas DataFrame with `smiles`, `score` columns, and compute
     the compound fingerprints when each batch is collated together

     a) The `__init__` function should take parameters

         * df: `pandas.DataFrame` object
         * smiles_col: string name of the column containing the compound smiles
         * score_col: string name of the column containing the score
         * fingerprint_fn: function that converts a list of ingerprints to an `numpy.array`

     b) The `__getitem__` function that given an index returns a `(smiles, score)` tuple

     c) The `collate_fn` that that takes a list of `(smiles, score)` tuples and produces
        a `(fingerprints, scores)` tuple, where fingerprints and scores are of type `torch.tensor`.
        Here you can use e.g. MolFeat or RDkit to do the featurization.

      d) Below we will use the FingerprintDataset when we create dataloader like this

        train_dataset = FingerprintDataset(
            df = train_df,
            smiles_col = smiles_col,
            score_col = score_col)
        train_loader = DataLoader(
            train_dataset,
            batch_size = batch_size,
            shuffle = True,
            collate_fn = train_dataset.collate_fn)

  3) Implement a step function that takes in the following parameters

       * model: a `torch.nn.Module` object
       * data_loader: a `torch.utils.data.DataLoader` object
       * optimizer: a `torch.optim.Optimizer` object
       * criterion: a function that takes in predictions and scores and produces a loss
       * device: a string name for the device
       * train: a boolean value indicating if run is in train or eval mode

     The function should iterate through all the batches of the `data_loader`, in each iteration

       * move the batch of fingerprints and scores to the GPU device
       * if in training mode set the optimizer gradient to zero
       * evaluate model on the fingerprints tensor to generate predicted scores
       * compute the loss comparing the predicted scores vs the given scores using the criterion
       * if in training mode, compute gradients on the loss as `loss.backwards()` and take a
         step of the optimizer as `optimizer.step()`.
       * returns the loss

  4) Implement a main function that does the following

     a) parses command line arguments

        * `train_path`, `val_path`, `test_path`: paths for to `.tsv` or `.parquet` files
        * `smiles_col`, `score_col`: column names of the smiles and score in the data files
        * `fingerprint_type`: a MolFeat fingerprint type
        * `hidden_dim`, `n_layers`: hyper-parameters of the `FingerprintNN`
        * `batch_size`, `learning_rate`, `n_epochs`: hyperparameters of model training
        * `wandb_project`: Weights and Biases project ID
        * `wandb_tags`: a list of tags to label the training run
        * device: the name of the compute device to use (e.g. 'cuda')

     B) Set up for model training

        * Initialize Weights and Biases for the project and given tags and pass in
          the command line arguments as the config with `config=vars(args)`.
        * load the training, validation, and test data from the provided paths. For each,
          create a `torch.utils.data.DataLoader` by passing in an instance of the `FingerprintDataset`
          and specifying the `batch_size`, `suffle` and `collate_fn` arguments.
        * Initialize the model and optimizer, moving the model to the device. You can use the `Adam`
          optimizer with `model.parameters()` the supplied liearning rate.

     C) Iterate over the requested number of epochs of the training set (once through the dataset),

        * Iterate over the batches in from the dataloader for the training set
          and interate over the batches from the validation dataloader and log the average loss for both
          to weights and biases

     D) At the end iterate over the test loss and log the average loss to weights and biases


### 3 Check and Refine the model till the is training
Getting the a new model to train is a bit of an art.

   1) Train the model on a dataset consisting of a single data point and check that it gives no errors
      and the loss goes to zero

   2) Add [dropout](https://www.jmlr.org/papers/v15/srivastava14a.html) to make the model more robust

      * Add a float `dropout` argument to `FingerprintNN`, specifying the fraction of weights each
        layer to mask during training. Have this argument passed in from the command line
      * After each `ReLU` layer, add a `torch.nn.Dropout(dropout)` layer
      * In the training loop, before iterating through the training dataset call `model.train()`
        and then before iterating through the validation or test datasets call `model.eval()`,
        this will tell the model to use drop out or not.

   3) Set up sensible baseline model and training hyper parameters on small representive dataset

      a) Sample data down to `1_000` each for the training, validation, and test datasets, and use
         ECFP fingerprints from MolFeat..

      b) Test baseline model parmaeters

          * `hidden_dim=256`
          * `n_layers=3`
          * `drop_out=0.1`
          * `batch_size=512`
          * `learning_rate=1e-3`
          * `epochs=100`

      c) Training requires getting the right amount of stochasticity to be able to escape local minima
         but not have exploding gradients. The `hidden_dims`, initialization strategy, `batch_size`
         and `learning_rate` all contribute to the stochasticity. A systematic way to set these is to
         given a network and initialization strategy, adjust the batch size to be as big as possible to
         fit into the gpu memory (which lowers the stochasticity) and then increase the learning rate
         (which raises the stochasticity) until the gradients explode, then back off a bit.

      d) To visually evaluate if the model is training, in Weights and Biases, look at the
         the training and validation losses progress over training. They should fall off in
         exponetial decay like curves. If something else is happening this can indicate problems
         that can be investigated further to improve model training.


### 4 Investigate the effect of hyper-parameter optimization
Hyper parameter optimization can affect not only the final predictive accuracy,
but the shape and behavior of training.

   1) Sweep each of the key hyperparameters, but feel free to go larger or smaller

      * hidden_dim: `[64, 128, 512, 1024, 2048]`
      * n_layers: `[1, 3, 5, 7]`
      * batch_size: `[64, 128, 512, 1024, 2048]`
      * learning_rate: `[1, 3, 10, 30, 100, 300] * 1e-4`

   2) Combine the best results from each sweep and then scale the dataset size to `100_000`
      ligands.

   3) Use the Weights and Biases API to download an plot the final train/valid/test losses of runs
      for each of the sweeps you ran.

          import wandb
          import pandas as pd

          api = wandb.Api()
          runs = api.runs("maomlab/BIOINF595w26_lab_08")

          all_summaries = []
          for run in runs:
              run_data = run.config | dict(run.summary)
              run_data['project'] = run.project
              run_data['name'] = run.name
              run_data['id'] = run.id
              all_summaries.append(run_data)
      
          all_summaries = pd.DataFrame(all_summaries)

### Questions

   a) Did you see different behaviors in the loss curves as across different parameters?
   b) Which hyperparameters had larger or smaller effects?
   c) What additional improvements/"tricks" would be good to try next to improve the quality of the model?

Please turn

   * A plot from Weights and Biases of representitive loss curves
   * The short answer + plots of the hyperparameter optimization
   * The code for `model.py` and hyperparameter optimization
