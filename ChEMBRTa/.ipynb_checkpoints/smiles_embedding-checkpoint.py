from metaflow import FlowSpec, Parameter, step, pypi_base

@pypi_base(packages={"transformers":"4.48.2",
                    "torch":"2.6.0",
                    "pandas":"2.2.3",
                    "numpy":"1.26.4",
                    }, python="3.10.13")
class SMILESFeatures(FlowSpec):
    file = Parameter("input", help="SMILES file")
    model_name = Parameter("model", help= "Model name", type=str)
    num_bacthes = Parameter("num_batch", help = "Number of SMILES data batches/chunks", type=int)
    out_dir = Parameter("output", help = "Output dir")

    @step
    def start(self):
        import os
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        print("Init SMILES representation learning")
        self.next(self.load_model, self.read_smiles)
    
    
    @step
    def load_model(self):
        from transformers import AutoModel, AutoTokenizer
        import torch 

        if torch.cuda.is_available():
            self.device = "cuda:0"
        else:
            self.device = "cpu"
        print("Working on: ",self.device)
        
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        print(type(self.tokenizer))

        self.model = AutoModel.from_pretrained(self.model_name)
        self.model = self.model.to(torch.device(self.device))
        self.next(self.join)

    @step
    def read_smiles(self):
        import pandas as pd
        self.smiles = pd.read_csv(self.file)
        self.smiles = self.smiles.dropna() #There are None SMILES
        self.next(self.join)
    
    @step
    def join(self, inputs):
        #Not possible to use varialbes in steps that are not the previous. So, flow is split, and model, tokenizer, and data are join together for further steps
        self.model = inputs.load_model.model
        self.tokenizer = inputs.load_model.tokenizer
        self.data = inputs.read_smiles.smiles
        self.device = inputs.load_model.device
        self.next(self.split_smiles)

    @step
    def split_smiles(self):
        import numpy as np

        self.data_arrays = np.array_split(self.data["norm_smiles"], self.num_bacthes)
        self.index_data_arrays = list(enumerate(self.data_arrays))
        self.next(self.gen_embeddings, foreach = "index_data_arrays", )

    @step
    def gen_embeddings(self):
        import torch
        index, data_batch = self.input
        self.batch_token = self.tokenizer(data_batch.to_list(),
                                          return_tensors = "pt",
                                          padding="max_length",              # Add special padding tokens to shorter sequences, ensuring fixed-size tensors 
                                          truncation=True,
                                          max_length=512 ).to(torch.device(self.device))
        with torch.no_grad():
            self.batch_token_embeddings = self.model(**self.batch_token)
        self.batch_token_embeddings_lhs = self.batch_token_embeddings.last_hidden_state.cpu()
        torch.save(self.batch_token_embeddings_lhs[:,0,:].clone(), self.out_dir + "/cls_smiles_batch_{0}".format(index))
        #torch.save(torch.mean(self.batch_token_embeddings_lhs, dim = 1), self.out_dir + "/mean_smiles_batch_{0}".format(index))

        del self.batch_token_embeddings_lhs
        del self.batch_token_embeddings
        del self.batch_token

        self.next(self.free_cuda)
    
    @step
    def free_cuda(self, _):
        import torch

        torch.cuda.empty_cache()
        
        self.next(self.end)


    @step
    def end(self):
        pass
    
if __name__ == "__main__":
    SMILESFeatures()
