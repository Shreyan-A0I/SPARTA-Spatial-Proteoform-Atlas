#!/bin/bash
# Script to run the SPARTA dashboard

echo "Starting SPARTA Dashboard..."
echo "The dashboard will open in your browser automatically."
echo "Press Ctrl+C to stop the server."
echo ""

cd "$(dirname "$0")"
conda run -n sparta streamlit run src/frontend/dashboard.py

