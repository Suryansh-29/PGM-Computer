import tkinter as tk
from tkinter import messagebox, filedialog, scrolledtext
import numpy as np

class PGMComputer:
    def __init__(self, root):
        self.root = root
        self.root.title("PGM Computer")
        self.root.geometry("750x700")
        self.root.configure(bg="#e0f7fa")  # light blue background
        
        # Frame for prime input
        top_frame = tk.Frame(root, bg="#e0f7fa")
        top_frame.pack(pady=10)
        tk.Label(top_frame, text="Enter prime number q:", bg="#e0f7fa", font=("Arial", 12, "bold")).grid(row=0, column=0, padx=5)
        self.q_entry = tk.Entry(top_frame, width=10, font=("Arial", 12))
        self.q_entry.grid(row=0, column=1, padx=5)
        self.q_entry.insert(0, "3")  # default q = 3
        
        # Frame for codewords input
        input_frame = tk.Frame(root, bg="#e0f7fa")
        input_frame.pack(pady=10)
        tk.Label(input_frame, text="Enter codewords (one per line, space separated):", 
                 bg="#e0f7fa", font=("Arial", 12, "bold")).pack(anchor="w", padx=5)
        self.input_text = scrolledtext.ScrolledText(input_frame, width=70, height=10, 
                                                     font=("Consolas", 12), bg="#ffffff", fg="#000000")
        self.input_text.pack(padx=5, pady=5)
        # Default example: a basis that is not already in standard form
        # This basis (over GF(3)) generates a [4,2]_3 code with minimum Hamming distance = 2.
        example = "1 0 1 0\n0 1 2 1\n1 1 1 2"
        self.input_text.insert(tk.END, example)
        
        # Frame for buttons
        btn_frame = tk.Frame(root, bg="#e0f7fa")
        btn_frame.pack(pady=10)
        self.compute_button = tk.Button(btn_frame, text="Compute PGM", command=self.compute_pgm, 
                                        font=("Arial", 12), bg="#00796b", fg="white")
        self.compute_button.grid(row=0, column=0, padx=10)
        self.clear_button = tk.Button(btn_frame, text="Clear Input", command=self.clear_input, 
                                      font=("Arial", 12), bg="#c62828", fg="white")
        self.clear_button.grid(row=0, column=1, padx=10)
        self.save_button = tk.Button(btn_frame, text="Save Output", command=self.save_output, 
                                     font=("Arial", 12), bg="#1565c0", fg="white")
        self.save_button.grid(row=0, column=2, padx=10)
        
        # Output text area
        output_frame = tk.Frame(root, bg="#e0f7fa")
        output_frame.pack(pady=10)
        tk.Label(output_frame, text="Output:", bg="#e0f7fa", font=("Arial", 12, "bold")).pack(anchor="w", padx=5)
        self.output_text = scrolledtext.ScrolledText(output_frame, width=70, height=20, 
                                                      font=("Consolas", 12), bg="#f1f8e9", fg="#1b5e20")
        self.output_text.pack(padx=5, pady=5)
    
    def clear_input(self):
        self.input_text.delete("1.0", tk.END)
    
    def get_input_codewords(self):
        data = self.input_text.get("1.0", tk.END).strip()
        if not data:
            return []
        codewords = []
        for line in data.splitlines():
            if line.strip():
                try:
                    # Convert each space-separated element to int
                    codewords.append(list(map(int, line.strip().split())))
                except Exception as e:
                    messagebox.showerror("Error", f"Invalid codeword line: '{line}'.\n{e}")
                    return []
        return codewords

    def compute_pgm(self):
        # Get q from entry
        try:
            q = int(self.q_entry.get())
            if q < 2:
                raise ValueError("q must be at least 2")
        except Exception as e:
            messagebox.showerror("Error", f"Invalid prime number q: {e}")
            return

        if not self.is_prime(q):
            messagebox.showerror("Error", f"{q} is not a prime number!")
            return

        # Get codewords from text input
        self.codewords = self.get_input_codewords()
        if not self.codewords:
            messagebox.showwarning("Warning", "No valid codewords input!")
            return

        # Work in GF(q)
        codewords_matrix = np.array(self.codewords) % q
        n = codewords_matrix.shape[1]
        # Note: np.linalg.matrix_rank is not GF(q)-aware, but if the basis is chosen properly it works for small examples.
        k = np.linalg.matrix_rank(codewords_matrix)
        if k == 0:
            messagebox.showerror("Error", "Invalid codewords! Rank is zero.")
            return

        # Compute systematic generator matrix in GF(q)
        G_sys, pivot_cols, non_pivot_cols = self.get_systematic_generator_matrix(codewords_matrix, q)
        if G_sys is None:
            messagebox.showerror("Error", "Could not compute systematic generator matrix.")
            return

        # The systematic generator matrix (PGM) has dimension k x n
        PGM = G_sys % q
        
        # Compute message-codeword mapping using the generator matrix
        mapping = self.get_message_codeword_mapping(PGM, q)
        min_distance = self.get_minimum_distance(mapping)

        output = f"Prime number q: {q}\n"
        output += f"Original Codewords (mod {q}):\n{codewords_matrix}\n\n"
        output += f"Systematic Generator Matrix (PGM) [k x n]:\n{PGM}\n\n"
        output += "Message-Codeword Mapping:\n"
        for msg, cw in mapping.items():
            output += f"Message {msg} -> Codeword {cw}\n"
        output += f"\nMinimum Hamming Distance: {min_distance}\n"

        self.output_text.delete("1.0", tk.END)
        self.output_text.insert(tk.END, output)

    def is_prime(self, num):
        if num < 2:
            return False
        for i in range(2, int(num**0.5) + 1):
            if num % i == 0:
                return False
        return True

    def mod_inverse(self, a, q):
        a = a % q
        for x in range(1, q):
            if (a * x) % q == 1:
                return x
        raise ValueError("Modular inverse does not exist.")

    def gauss_elim_modq(self, matrix, q):
        A = matrix.copy() % q
        rows, cols = A.shape
        pivot_row = 0
        pivot_cols = []
        for col in range(cols):
            if pivot_row >= rows:
                break
            pivot_found = False
            for r in range(pivot_row, rows):
                if A[r, col] % q != 0:
                    if r != pivot_row:
                        A[[pivot_row, r]] = A[[r, pivot_row]]
                    pivot_found = True
                    break
            if not pivot_found:
                continue
            pivot = A[pivot_row, col] % q
            inv = self.mod_inverse(pivot, q)
            A[pivot_row] = (A[pivot_row] * inv) % q
            pivot_cols.append(col)
            for r in range(rows):
                if r != pivot_row:
                    factor = A[r, col] % q
                    A[r] = (A[r] - factor * A[pivot_row]) % q
            pivot_row += 1
        return A, pivot_cols

    def get_systematic_generator_matrix(self, matrix, q):
        R, pivot_cols = self.gauss_elim_modq(matrix, q)
        k = len(pivot_cols)
        if k == 0:
            return None, None, None
        rows, cols = R.shape
        G_temp = R[:k, :].copy()
        non_pivot_cols = [i for i in range(cols) if i not in pivot_cols]
        permutation = pivot_cols + non_pivot_cols
        G_sys = G_temp[:, permutation]
        I_k = np.eye(k, dtype=int) % q
        # Make sure the left block is the identity matrix by eliminating extra entries.
        if not np.array_equal(G_sys[:, :k] % q, I_k):
            for i in range(k):
                if G_sys[i, i] % q == 0:
                    for j in range(i+1, k):
                        if G_sys[j, i] % q != 0:
                            G_sys[[i, j]] = G_sys[[j, i]]
                            break
                for j in range(k):
                    if j != i and G_sys[j, i] % q != 0:
                        factor = G_sys[j, i] % q
                        inv = self.mod_inverse(G_sys[i, i] % q, q)
                        G_sys[j] = (G_sys[j] - factor * inv * G_sys[i]) % q
        return G_sys % q, pivot_cols, non_pivot_cols

    def get_message_codeword_mapping(self, G_sys, q):
        k = G_sys.shape[0]
        n = G_sys.shape[1]
        mapping = {}
        for i in range(q**k):
            msg = self.int_to_base_array(i, k, q)
            codeword = np.mod(np.dot(msg, G_sys), q)
            mapping[''.join(map(str, msg))] = ''.join(map(str, codeword))
        return mapping

    def int_to_base_array(self, number, length, base):
        arr = [0] * length
        for i in range(length - 1, -1, -1):
            arr[i] = number % base
            number //= base
        return np.array(arr, dtype=int)

    def get_minimum_distance(self, mapping):
        distances = []
        keys = list(mapping.keys())
        # Exclude the zero message (e.g., "000") when computing minimum nonzero weight.
        nonzero_keys = [key for key in keys if any(ch != '0' for ch in key)]
        for i in range(len(nonzero_keys)):
            cw1 = np.array(list(map(int, mapping[nonzero_keys[i]])))
            for j in range(i+1, len(nonzero_keys)):
                cw2 = np.array(list(map(int, mapping[nonzero_keys[j]])))
                distances.append(np.sum(cw1 != cw2))
        return min(distances) if distances else 0

    def save_output(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".txt",
                                                 filetypes=[("Text Files", "*.txt")])
        if not file_path:
            return
        try:
            with open(file_path, 'w') as file:
                file.write(self.output_text.get("1.0", tk.END))
            messagebox.showinfo("Success", "Output saved successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save file: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = PGMComputer(root)
    root.mainloop()
