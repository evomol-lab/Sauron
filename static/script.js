document.addEventListener('DOMContentLoaded', () => {
    const themeToggle = document.getElementById('theme-toggle');
    const logoImg = document.getElementById('logo-img');
    const themeIconLight = document.getElementById('theme-icon-light');
    const themeIconDark = document.getElementById('theme-icon-dark');
    
    // Check saved theme
    if (localStorage.getItem('theme') === 'light') {
        document.body.classList.add('light-mode');
        logoImg.src = '/static/sauron-rins-logo.png';
        themeIconDark.style.display = 'none';
        themeIconLight.style.display = 'block';
    }

    themeToggle.addEventListener('click', () => {
        document.body.classList.toggle('light-mode');
        if (document.body.classList.contains('light-mode')) {
            localStorage.setItem('theme', 'light');
            logoImg.src = '/static/sauron-rins-logo.png';
            themeIconDark.style.display = 'none';
            themeIconLight.style.display = 'block';
        } else {
            localStorage.setItem('theme', 'dark');
            logoImg.src = '/static/sauron-rins-logo-inv.png';
            themeIconLight.style.display = 'none';
            themeIconDark.style.display = 'block';
        }
    });

    const dropZone = document.getElementById('drop-zone');
    const fileInput = document.getElementById('file-input');
    const fileNameDisplay = document.getElementById('file-name');
    const submitBtn = document.getElementById('submit-btn');
    const uploadForm = document.getElementById('upload-form');
    const errorMessage = document.getElementById('error-message');
    const successState = document.getElementById('success-state');
    const downloadBtn = document.getElementById('download-btn');
    const resetBtn = document.getElementById('reset-btn');

    let selectedFile = null;

    // Handle Drag & Drop
    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
        dropZone.addEventListener(eventName, preventDefaults, false);
    });

    function preventDefaults(e) {
        e.preventDefault();
        e.stopPropagation();
    }

    ['dragenter', 'dragover'].forEach(eventName => {
        dropZone.addEventListener(eventName, () => dropZone.classList.add('dragover'), false);
    });

    ['dragleave', 'drop'].forEach(eventName => {
        dropZone.addEventListener(eventName, () => dropZone.classList.remove('dragover'), false);
    });

    dropZone.addEventListener('drop', (e) => {
        const dt = e.dataTransfer;
        const files = dt.files;
        handleFiles(files);
    });

    // Handle Click Upload
    dropZone.addEventListener('click', () => fileInput.click());
    
    fileInput.addEventListener('change', function() {
        handleFiles(this.files);
    });

    function handleFiles(files) {
        if (files.length > 0) {
            const file = files[0];
            const validExtensions = ['.pdb', '.cif'];
            const fileExt = file.name.substring(file.name.lastIndexOf('.')).toLowerCase();
            
            if (validExtensions.includes(fileExt)) {
                selectedFile = file;
                fileNameDisplay.textContent = file.name;
                submitBtn.disabled = false;
                errorMessage.style.display = 'none';
            } else {
                fileNameDisplay.textContent = 'Invalid file format. Please upload .pdb or .cif';
                selectedFile = null;
                submitBtn.disabled = true;
            }
        }
    }

    // Form Submission via Fetch API
    uploadForm.addEventListener('submit', async (e) => {
        e.preventDefault();
        
        if (!selectedFile) return;

        // UI State: Loading
        submitBtn.classList.add('loading');
        submitBtn.disabled = true;
        errorMessage.style.display = 'none';

        const formData = new FormData();
        formData.append('file', selectedFile);
        formData.append('add_h', document.getElementById('add_h').checked);
        formData.append('strict_angle', document.getElementById('strict_angle').checked);
        formData.append('remove_multiples', document.getElementById('remove_multiples').checked);

        try {
            const response = await fetch('/upload', {
                method: 'POST',
                body: formData
            });

            const data = await response.json();

            if (!response.ok) {
                throw new Error(data.error || 'Computation failed.');
            }

            // UI State: Success
            uploadForm.style.display = 'none';
            successState.style.display = 'block';
            downloadBtn.href = data.downloadUrl;

        } catch (error) {
            // UI State: Error
            errorMessage.textContent = error.message;
            errorMessage.style.display = 'block';
        } finally {
            submitBtn.classList.remove('loading');
            submitBtn.disabled = false;
        }
    });

    // Reset Form
    resetBtn.addEventListener('click', () => {
        successState.style.display = 'none';
        uploadForm.style.display = 'block';
        uploadForm.reset();
        selectedFile = null;
        fileNameDisplay.textContent = 'No file selected';
        submitBtn.disabled = true;
        errorMessage.style.display = 'none';
    });
});
