document.addEventListener("DOMContentLoaded", function() {
    document.querySelectorAll(".options-expander").forEach(expander => {
        expander.addEventListener("click", function() {
            expander.parentElement.querySelector(".options").classList.toggle("hide");
            expander.querySelector(".fa").classList.toggle("fa-chevron-down");
            expander.querySelector(".fa").classList.toggle("fa-chevron-up");
        });
    });

    document.querySelectorAll(".col-3 label").forEach(label => {
        label.addEventListener("click", function() {
            label.parentElement.querySelector("input").click();
        });
    });

    document.querySelectorAll("h2:not(.setting-group-title)").forEach(h2 => {
        let my_headings = document.querySelectorAll("h2.setting-group-title");
        for (let heading of my_headings) {
            if (heading.innerText == h2.innerText) {
                h2.style = "height: 0; margin: 0; overflow: hidden;";
                break;
            }
        }
    });
});