import tensorflow as tf
import keras
from keras import layers
import numpy as np

def forward_net():
    # input: sequence
    inputs = keras.Input(shape=(512,), dtype="int32") 
    x = layers.Embedding(
        input_dim=21, 
        output_dim=256, # output-dim of layer
        mask_zero=True
        )(inputs)
    x = layers.Dense(256, activation="relu")(x)
    outputs = layers.Dense(512, activation="linear")(x)

    model = keras.Model(inputs=inputs, outputs=outputs)
    return model

def run_model(model, train_list, val_list, test_list):
    model.compile (
        optimizer="rmsprop",
        loss="mse",
        metrics =["mae"]
        )
    
    train_samples, train_labels = zip(*train_list)
    train_samples = np.array(train_samples)
    train_labels = np.array(train_labels)

    val_samples, val_labels = zip(*val_list)
    val_samples = np.array(val_samples)
    val_labels = np.array(val_labels)

    test_samples, test_labels = zip(*test_list)
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
    

