window.addEventListener('load', function() {
    // Create a container for logos
    var logoDiv = document.createElement('div');
    logoDiv.className = 'multiple-logos';
    
    // Create logo images
    var logo1 = document.createElement('img');
    logo1.src = window.quartoLogos.logo1;
    logo1.alt = 'SPC Logo';
    
    var logo2 = document.createElement('img');
    logo2.src = window.quartoLogos.logo2;
    logo2.alt = 'NOAA Fisheries Logo';
    
    // Add logos to container
    logoDiv.appendChild(logo1);
    logoDiv.appendChild(logo2);
    
    // Add the logo container to the reveal container
    document.querySelector('.reveal').appendChild(logoDiv);
    
    // Update logo positioning on slide change
    function updateLogoPosition() {
        // This is handled by CSS classes now
        // Just ensuring logos are visible
        logoDiv.style.display = 'flex';
    }
    
    // Set up event listeners for slide changes
    Reveal.on('slidechanged', updateLogoPosition);
    
    // Initial positioning
    Reveal.on('ready', updateLogoPosition);
});
