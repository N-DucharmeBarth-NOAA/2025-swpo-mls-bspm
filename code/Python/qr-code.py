import qrcode
from PIL import Image
import os
import sys

def generate_qr_code(url, filename=None, size=10, border=4):
    """
    Generate a QR code for a given URL
    
    Args:
        url (str): The website URL to encode
        filename (str): Output filename (optional, defaults to 'qr_code.png')
        size (int): Size of each box in pixels (default: 10)
        border (int): Border size in boxes (default: 4)
    """
    
    # Add https:// if no protocol is specified
    if not url.startswith(('http://', 'https://')):
        url = 'https://' + url
    
    # Create QR code instance
    qr = qrcode.QRCode(
        version=1,  # Controls size (1 is smallest)
        error_correction=qrcode.constants.ERROR_CORRECT_L,  # About 7% error correction
        box_size=size,
        border=border,
    )
    
    # Add data and optimize
    qr.add_data(url)
    qr.make(fit=True)
    
    # Create image
    img = qr.make_image(fill_color="black", back_color="white")
    
    # Set filename if not provided
    if filename is None:
        # Create safe filename from URL
        safe_url = url.replace('https://', '').replace('http://', '').replace('/', '_').replace(':', '')
        filename = f"qr_{safe_url}.png"
    
    # Save image
    img.save(filename)
    print(f"QR code saved as: {filename}")
    print(f"QR code contains: {url}")
    
    return filename

def batch_generate_qr_codes(urls, output_dir="qr_codes"):
    """
    Generate QR codes for multiple URLs
    
    Args:
        urls (list): List of URLs
        output_dir (str): Directory to save QR codes
    """
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    generated_files = []
    
    for i, url in enumerate(urls):
        # Create safe filename
        safe_url = url.replace('https://', '').replace('http://', '').replace('/', '_').replace(':', '')
        filename = os.path.join(output_dir, f"qr_{i+1}_{safe_url}.png")
        
        try:
            generate_qr_code(url, filename)
            generated_files.append(filename)
        except Exception as e:
            print(f"Error generating QR code for {url}: {e}")
    
    return generated_files

def main():
    """Interactive command-line interface"""
    
    print("QR Code Generator for Websites")
    print("=" * 35)
    
    while True:
        print("\nOptions:")
        print("1. Generate single QR code")
        print("2. Generate multiple QR codes")
        print("3. Exit")
        
        choice = input("\nEnter your choice (1-3): ").strip()
        
        if choice == "1":
            url = input("Enter website URL: ").strip()
            if url:
                try:
                    filename = input("Enter filename (press Enter for default): ").strip()
                    if not filename:
                        filename = None
                    
                    generate_qr_code(url, filename)
                except Exception as e:
                    print(f"Error: {e}")
            else:
                print("Please enter a valid URL.")
        
        elif choice == "2":
            print("Enter URLs (one per line, empty line to finish):")
            urls = []
            while True:
                url = input().strip()
                if not url:
                    break
                urls.append(url)
            
            if urls:
                try:
                    output_dir = input("Enter output directory (press Enter for 'qr_codes'): ").strip()
                    if not output_dir:
                        output_dir = "qr_codes"
                    
                    generated_files = batch_generate_qr_codes(urls, output_dir)
                    print(f"\nGenerated {len(generated_files)} QR codes in '{output_dir}' directory")
                except Exception as e:
                    print(f"Error: {e}")
            else:
                print("No URLs provided.")
        
        elif choice == "3":
            print("Goodbye!")
            break
        
        else:
            print("Invalid choice. Please enter 1, 2, or 3.")

if __name__ == "__main__":
    # Check if URL is provided as command line argument
    if len(sys.argv) > 1:
        url = sys.argv[1]
        filename = sys.argv[2] if len(sys.argv) > 2 else None
        generate_qr_code(url, filename)
    else:
        main()

# Example usage:
# python qr_generator.py https://www.google.com
# python qr_generator.py github.com my_qr.png