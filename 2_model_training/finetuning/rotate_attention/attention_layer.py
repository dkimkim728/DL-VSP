from torch.nn import Linear, Module
from fast_transformers.attention import AttentionLayer
from fast_transformers.events import EventDispatcher, QKVEvent
from .rotary import RotaryEmbedding, apply_rotary_pos_emb
from fast_transformers.masking import LengthMask
import torch

class RotateAttentionLayer(AttentionLayer):
    """Rotate attention layer inherits from fast_transformer attention layer. 
        The only thing added is an Embedding encoding, for more information
        on the attention layer see the fast_transformers code.
    """
    def __init__(self, attention, d_model, n_heads, d_keys=None, d_values=None, event_dispatcher=""):
        super(RotateAttentionLayer, self).__init__(attention, d_model, n_heads, d_keys=d_keys, d_values=d_values, event_dispatcher=event_dispatcher)
        self.rotaryemb = RotaryEmbedding(d_keys)
        print('Using Rotation Embedding')

    def forward(self, queries, keys, values, attn_mask, query_lengths, key_lengths):
        """
        Forward pass through the rotate attention layer. It applies rotary embeddings to the queries and keys.
        """
        # Extract the dimensions into local variables
        if isinstance(key_lengths, LengthMask):
            key_lengths_tensor = key_lengths.lengths  # Use the 'lengths' attribute
            key_lengths = torch.tensor(key_lengths_tensor).to(torch.int64)

            # Add the print statements and the key-length check here
            print("keys shape:", keys.shape)  # Print the shape of keys
            print("key_lengths shape:", key_lengths.shape)  # Print the shape of key_lengths

            if keys.size(1) < key_lengths_tensor.max().item():  # Compare with max key length
                raise ValueError("Key length mismatch")
        
        # Proceed with the original logic
        N, L, _ = queries.shape
        _, S, _ = keys.shape
        H = self.n_heads

        # Project the queries/keys/values
        queries = self.query_projection(queries).view(N, L, H, -1)
        keys = self.key_projection(keys).view(N, S, H, -1)
        cos, sin = self.rotaryemb(queries)
        queries, keys = apply_rotary_pos_emb(queries, keys, cos, sin)
        values = self.value_projection(values).view(N, S, H, -1)

        # Debugging: Print tensor shapes
        print(f"queries shape: {queries.shape}, keys shape: {keys.shape}, values shape: {values.shape}")
        print(f"key_lengths shape: {key_lengths.shape}")

        # Let the world know of the qkv
        self.event_dispatcher.dispatch(QKVEvent(self, queries, keys, values))

        # Compute the attention
        new_values = self.inner_attention(
            queries,
            keys,
            values,
            attn_mask,
            query_lengths,
            key_lengths  # Use the key_lengths tensor directly
        ).view(N, L, -1)

        # Project the output and return
        return self.out_projection(new_values)