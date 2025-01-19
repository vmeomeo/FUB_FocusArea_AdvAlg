import argparse


def is_prime(n):
    if n <= 1:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True


def compute_primes(length):
    primes = []
    num = 2
    while len(primes) < length:
        if is_prime(num):
            primes.append(num)
        num += 1
    return primes


def main():
    parser = argparse.ArgumentParser(description="Compute a specified number of prime numbers.")
    parser.add_argument("-pathtofile", type=str, required=True,
                        help="Path to the file containing the number of primes to compute")
    args = parser.parse_args()

    try:
        with open(args.pathtofile, "r") as file:
            length = int(file.read().strip())
            primes = compute_primes(length)
            print(f"The first {length} prime numbers are: {primes}")
    except FileNotFoundError:
        print(f"Error: The file '{args.pathtofile}' was not found.")
    except ValueError:
        print(f"Error: The file '{args.pathtofile}' does not contain a valid integer.")


if __name__ == "__main__":
    main()
