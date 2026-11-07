# Adapted from the user's original metadata_qc_helpers final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
haversine_km <- function(lat1, lon1, lat2, lon2) {
    earth_radius_km <- 6371.0088
    to_rad <- pi/180
    phi1 <- lat1 * to_rad
    phi2 <- lat2 * to_rad
    delta_phi <- (lat2 - lat1) * to_rad
    delta_lambda <- (lon2 - lon1) * to_rad
    a <- sin(delta_phi/2)^2 + cos(phi1) * cos(phi2) * sin(delta_lambda/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    earth_radius_km * c
}
