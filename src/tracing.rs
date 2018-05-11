#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Tracing {
    Bidirectional,
    Importance,
    Light,
}
