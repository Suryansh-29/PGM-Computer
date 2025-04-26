# PGM Computer

A GUI-based tool to compute the **Systematic Generator Matrix (PGM)** for a given set of codewords over a finite field **GF(q)**.
Developed using **Python**, **Tkinter**, and **NumPy**.

---

## ✨ Features

- Input codewords manually (with space-separated values).
- Select the prime number \( q \) (works over **GF(q)**).
- Computes:
  - **Systematic Generator Matrix** (PGM) using modular Gaussian elimination.
  - **Message to Codeword Mapping**.
  - **Minimum Hamming Distance** of the code.
- Save output to a `.txt` file.
- Clean, user-friendly Tkinter GUI.
- Built-in examples for easy testing.

---

## 💻 Technologies Used

- Python 3
- Tkinter (for GUI)
- NumPy (for matrix operations)

---


## 🚀 How to Run

1. Make sure you have **Python 3.x** installed.
2. Install NumPy (if not already installed):
   ```bash
   pip install numpy
   ```
3. Run the script:
   ```bash
   python pgm_computer.py
   ```

---

## 🧐 How It Works

- You input the codewords and prime number \( q \).
- The program performs **Gaussian elimination mod q** to find a systematic form.
- Displays the generator matrix, mapping between messages and codewords, and minimum Hamming distance.

---



## 📝 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## ✍️ Author

- **Suryansh Humane**
  - [GitHub](https://github.com/Suryansh-29)

---
