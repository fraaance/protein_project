import tensorflow as tf
import keras
from keras import layers
import numpy as np
import matplotlib.pyplot as plt

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
    x = layers.Dense(256, activation="relu")(x)
    outputs = layers.Dense(512, activation="linear")(x)

    model = keras.Model(inputs=inputs, outputs=outputs)
    return model

def run_model(model, train_list, val_list, test_list, distances_path):
    # update:
    # train_list, val_list, test_list only contain seq_names, aa_encoded

    model.compile(
        optimizer="rmsprop",
        loss=loss_mask_mse,
        metrics=[loss_mask_mae],
        run_eagerly=True
    )

    train_dataset = tf_set(train_list, distances_path)
    val_dataset = tf_set(val_list, distances_path)
    test_dataset = tf_set(test_list, distances_path)

    # run model: 
    history = model.fit(
        train_dataset,
        validation_data=val_dataset,
        epochs=10
    )

    test_loss, test_mae = model.evaluate(test_dataset)
    print(f" model trained with train set of size: {len(train_list)} ".center(num_stars, "*"))
    
    epochs = range(1, 11)
    print(f"test MAE: {test_mae}")

    fig, ax = plt.subplots()
    ax.plot(epochs, history.history["loss"], "bo-", label="test loss")
    ax.plot(epochs, history.history["val_loss"], "ro-", label="val loss")
    ax.set_title("Training Loss")
    ax.set_xlabel("Epochs")
    ax.set_ylabel("Loss")
    ax.legend()
    return model, fig


def tf_set(list, distances_path):
    dataset = tf.data.Dataset.from_generator(
        lambda: batch_generator(list, distances_path),
        output_signature=(
            tf.TensorSpec(shape=(512,), dtype=tf.int32),
            tf.TensorSpec(shape=(512, 512, 2),dtype=tf.float32)
        )
    )
    return dataset.batch(8).prefetch(1)


def batch_generator(data_list, distances_path):
    # returns the sequence as input, and (dist_matrix, mask) as output
    for seq_name, seq in data_list:
        file_path = distances_path / f"{seq_name}.npz"

        if file_path.exists():
            with np.load(file_path) as data:
                dist = data["dist_matrix"].astype(np.float32)
                mask = data["mask"].astype(np.float32)

            y = np.stack([dist, mask], axis=-1)
            yield seq, y #return seq, y

def loss_mask_mse(y_true, y_pred):
    dist = y_true[..., 0]
    mask = y_true[..., 1]

    sq_error = tf.square(dist - y_pred) * mask
    return tf.reduce_sum(sq_error) / tf.maximum(tf.reduce_sum(mask), 1.0)

def loss_mask_mae(y_true, y_pred):
    dist = y_true[..., 0]
    mask = y_true[..., 1]
    abs_error = tf.abs(dist - y_pred) * mask

    return tf.reduce_sum(abs_error) / tf.maximum(tf.reduce_sum(mask), 1.0)