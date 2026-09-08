import tensorflow as tf
import keras
from keras import layers
import numpy as np
num_stars = 70

def forward_net():
    # input: sequence
    inputs = keras.Input(shape=(512,), dtype="int32") 
    x = layers.Embedding(
        input_dim=21, 
        output_dim=256, # output-dim of layer
        #mask_zero=True
        )(inputs)
    x = layers.Dense(256, activation="relu")(x)
    outputs = layers.Dense(512, activation="linear")(x)

    model = keras.Model(inputs=inputs, outputs=outputs)
    return model

def run_model(model, train_list, val_list, test_list):
    model.compile (
        optimizer="rmsprop",
        loss=loss_mask_mse,
        metrics =[loss_mask_mae],
        run_eagerly=True
    )
    
    train_samples, train_labels = zip(*train_list)
    train_samples = np.array(train_samples)
    train_dist, train_mask = unzip_dist_matrix(train_labels)

    val_samples, val_labels = zip(*val_list)
    
    val_samples = np.array(val_samples)
    val_dist, val_mask = unzip_dist_matrix(val_labels)

    test_samples, test_labels = zip(*test_list)
    test_samples = np.array(test_samples)
    test_dist, test_mask = unzip_dist_matrix(test_labels)

    train_y = np.stack([train_dist, train_mask], axis=-1)
    val_y = np.stack([val_dist, val_mask], axis=-1)
    test_y = np.stack([test_dist, test_mask], axis=-1)

    pred = model(train_samples[:32])

    model.fit(
        train_samples, train_y,
        validation_data=(val_samples, val_y),
        epochs=20
    )

    test_loss, acc = model.evaluate(test_samples, test_y) 
    
    print(f" model trained with train set of size: {len(train_samples)} ".center(num_stars, "*"))
    print(f"test acc: {acc}")

def unzip_dist_matrix(label_list):
    dist_list, mask_list = zip(*label_list)
    dist = np.stack(dist_list, axis=0)
    mask = np.stack(mask_list, axis=0)
    return dist, mask
    
def combine_dist_mask(dist_matrix, mask):
    return np.stack([dist_matrix, mask], axis=-1)

def loss_mask_mse(y_true, y_pred):
    # based chatgpt... need to improve
    dist = y_true[..., 0]
    mask = y_true[..., 1]

    sq_error = tf.square(dist - y_pred) * mask

    return tf.reduce_sum(sq_error) / tf.maximum(tf.reduce_sum(mask), 1.0)

def loss_mask_mae(y_true, y_pred):

    dist = y_true[..., 0]
    mask = y_true[..., 1]

    abs_error = tf.abs(dist - y_pred) * mask

    return tf.reduce_sum(abs_error) / tf.maximum(tf.reduce_sum(mask), 1.0)
    #model.save("/Users/franzweisel/Downloads/project/base_model.keras")
    