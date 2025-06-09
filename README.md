# Model-Free Adaptive Iterative Learning Integral Terminal Sliding Mode Control of Exoskeleton Robots

This repository contains MATLAB code and simulations for our publication:

📄 **Title**: Model-Free Adaptive Iterative Learning Integral Terminal Sliding Mode Control of Exoskeleton Robots  
📰 **Journal**: Journal of Vibration and Control, Vol. 28(21–22), 2022  
🔗 [DOI: 10.1177/10775463211026031](https://doi.org/10.1177/10775463211026031)

---

## 🧠 Abstract

This study proposes a model-free adaptive iterative learning control strategy based on integral terminal sliding mode control (ITSMC) for exoskeleton robots. Without requiring an explicit system model, the method ensures high-precision tracking and finite-iteration convergence using only measured input/output data. Robustness is enhanced via a sliding surface that eliminates chattering while maintaining smooth control signals. Simulation results verify superior performance over existing methods.

---

## 🛠 Usage

To run the MATLAB code:

1. Clone this repository to your local machine.
2. Open the `.m` simulation scripts in MATLAB.
3. Customize controller parameters if needed.
4. Run the scripts to visualize the tracking performance and convergence behavior.

> ✅ No commercial solver like MOSEK is required. YALMIP-compatible solvers such as **SDPT3** or **SeDuMi** are sufficient.

---

## 📜 License and Contact Info

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details. You are welcome to customize the controller parameters and reuse this code for your control system applications.

If you have any questions or encounter issues, please feel free to contact me.

Enjoy exploring and using the code!

---

## 📚 Citation

If you found this repository useful in your research, please cite the original paper as:

```bibtex
@article{esmaeili2022model,
  title={Model-free adaptive iterative learning integral terminal sliding mode control of exoskeleton robots},
  author={Esmaeili, Babak and Madani, Seyedeh Sepideh and Salim, Mina and Baradarannia, Mahdi and Khanmohammadi, Sohrab},
  journal={Journal of Vibration and Control},
  volume={28},
  number={21-22},
  pages={3120--3139},
  year={2022},
  publisher={SAGE Publications},
  doi={10.1177/10775463211026031}
}
