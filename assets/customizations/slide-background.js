// Refactored clean implementation for slide backgrounds
document.addEventListener('DOMContentLoaded', function() {
  // Function to apply appropriate classes to slide content and backgrounds
  function updateSlideStyles() {
    // Get current slide
    const currentSlide = Reveal.getCurrentSlide();
    if (!currentSlide) return;
    
    // Handle background classes
    document.body.classList.remove('wide-sidebar-active', 'narrow-sidebar-active');
    document.body.classList.remove('wide-sidebar-bg', 'narrow-sidebar-bg');
    
    // Remove all horizontal lines first
    document.querySelectorAll('.horizontal-line').forEach(line => line.remove());
    
    // First remove all content positioning classes from ALL slides to prevent inheritance issues
    document.querySelectorAll('.slides section').forEach(slide => {
      // Remove content positioning classes
      slide.querySelectorAll('h1, h2, h3, p, ul, ol, .columns, [class*="column"]').forEach(el => {
        el.classList.remove('content-with-wide-sidebar', 'content-with-narrow-sidebar');
      });
      
      // Remove background classes
      slide.classList.remove('wide-sidebar-bg', 'narrow-sidebar-bg');
    });
    
    // Apply appropriate styling based on slide type
    if (currentSlide.classList.contains('left-bg-image')) {
      // Title slide with wide sidebar
      document.body.classList.add('wide-sidebar-bg');
      document.body.classList.add('wide-sidebar-active');
      
      // Position content with wide margin
      currentSlide.querySelectorAll('h1, h2, h3, p, ul, ol, .columns, [class*="column"]').forEach(el => {
        el.classList.add('content-with-wide-sidebar');
      });
      
    } else if (currentSlide.classList.contains('narrow-left-bg-image')) {
      // Regular slide with narrow sidebar
      document.body.classList.add('narrow-sidebar-bg');
      document.body.classList.add('narrow-sidebar-active');
      
      // Position content with slight margin
      currentSlide.querySelectorAll('h1, h2, h3, p, ul, ol, .columns, [class*="column"]').forEach(el => {
        el.classList.add('content-with-narrow-sidebar');
      });
      
      // Add horizontal line directly to this slide
      const horizontalLine = document.createElement('div');
      horizontalLine.className = 'horizontal-line';
      currentSlide.appendChild(horizontalLine);
    }
  }
  
  // Update styles when slides change
  Reveal.on('slidechanged', updateSlideStyles);
  
  // Initial setup
  Reveal.on('ready', updateSlideStyles);
  
  // Override layout to maintain our styling
  const originalLayout = Reveal.layout;
  Reveal.layout = function() {
    originalLayout.apply(this, arguments);
    updateSlideStyles();
  };
});
