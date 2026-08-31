import tensorflow as tf
import keras
from keras import layers
import numpy as np

class FeedForwardNet:

    def __init__(self, train_list, val_list, test_list):
        self.train_list = train_list
        self.test_list = test_list
        self.val_list = val_list

    def forward_net(self):
        
        # input: sequence
        inputs = keras.Input(shape=(256,), dtype="int32") 
        x = layers.Embedding(
            input_dim=21, 
            output_dim=128, # output-dim of layer
            mask_zero=True
            )(inputs)
        x = layers.Dense(128, activation="relu")(x)
        outputs = layers.Dense(256, activation="linear")(x)

        model = keras.Model(inputs=inputs, outputs=outputs)
        return model

    def run_model(self, model):
        model.compile (
            optimizer="rmsprop",
            loss="mse",
            metrics =["mae"]
            )
        
        train_samples, train_labels = zip(*self.train_list)
        train_samples = np.array(train_samples)
        train_labels = np.array(train_labels)

        val_samples, val_labels = zip(*self.val_list)
        val_samples = np.array(val_samples)
        val_labels = np.array(val_labels)

        test_samples, test_labels = zip(*self.test_list)
        test_samples = np.array(test_samples)
        test_labels = np.array(test_labels)

        model.fit(
            train_samples, 
            train_labels, 
            validation_data=(val_samples, val_labels),
            epochs=10, 
            batch_size=32
        )

        test_loss, test_acc = model.evaluate(test_samples, test_labels)
        print(f"test accuracy: {test_acc}")
        

