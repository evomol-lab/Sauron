FROM python:3.10-slim

# Install system dependencies (e.g., dssp)
RUN apt-get update && apt-get install -y \
    dssp \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/sauron

# Copy requirements first for cache efficiency
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source code
COPY . .

# Create uploads directory (can be overridden by SAURON_UPLOAD_DIR env var)
RUN mkdir -p uploads && chmod 777 uploads

# Expose Flask port
EXPOSE 5000

# Run the Flask app
CMD ["python", "SauronGUI.py"]
