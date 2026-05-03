from quantum_tunnelling.simulation import run_demo


if __name__ == "__main__":
    output_dir = run_demo()
    print(f"Saved plots to {output_dir.resolve()}")
